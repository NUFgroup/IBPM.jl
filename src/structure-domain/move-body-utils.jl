"""
Utilities to support motion of a body. Three categories of utilities here:

    (3) Driver functions that move the body based on the motion type
        prescribed by the user
    (2) Helper functions and structs that convert user-prescribed info including
        rotation to motion at the discrete surface points needed by the code
    (1) Low-level functions for determining motion type, number of body points,
        and other stuff useful for the main functions
"""

"
Type (1): low-level functions
"
function MotionType( body::V where V<:Body )
    return typeof(body.motion)
end

function MotionType( bodies::Array{V, 1} where V<:Body )
    motions = [typeof(bodies[i].motion) for i=1:length(bodies)]
    if all(motions .== motions[1])
        return motions[1]
    else
        return Motion
    end
end

function BodyType( body::V where V<:Body )
    return typeof(body)
end

function BodyType( bodies::Array{V, 1} where V<:Body )
    btypes = [typeof(bodies[i]) for i=1:length(bodies)]
    if all(btypes .== btypes[1])
        return btypes[1]
    else
        return Body
    end
end

function get_body_info( bodies::Array{V, 1} where V <: Body )
    # determine the num of body points per body
    nf = [length(bodies[j].xb) for j=1:length(bodies)]
    nb = [size(bodies[j].xb, 1) for j=1:length(bodies)]
    return nb, nf
end

function get_ub(bodies::Array{V, 1} where V<:Body)
    nb, nf = get_body_info(bodies)
    ub = zeros(Float64, sum(nb), 2)

    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j=1:length(bodies)
        #move_body!(bodies[j].xb, bodies[j].ub, bodies[j].motion, t)  # Not time-dependent
        ub[nbod_tally .+ (1:nb[j]), :] .= bodies[j].ub

        # update body index
        nbod_tally += nb[j];
    end
    return [ub[:, 1]; ub[:, 2]]  # Stack [ub; vb]
end


"
Type (2): facilitate transfer of desired rotation motion to
    motion of surface points
"

#Struct to hold maps for position and velocity, modeled after Rowley's C++ code
#    position transformation:
#     x = pos .+ Rx*x
# velocity transformation:
#     v = vel .+ Rv*v

struct TangentSE2
    pos::Array{Float64, 1}   # [x, y] position of center
    vel::Array{Float64, 1}   # [ẋ, ẏ] velocity of center
    Rx::Any  # maps positions [x, y]
    Ru::Any  # maps velocities [ẋ, ẏ]
end


function get_transformation(
        motion::MotionFunction,
        t::Float64   # Should be time of the *next* step (i.e. t[k+1])
        )
    # x, y, θ = motion.xc(t)
    # ẋ, ẏ, θ̇ = motion.uc(t)

    x = motion.xc[1](t)
    y = motion.xc[2](t)
    θ = motion.xc[3](t)

    ẋ = motion.uc[1](t)
    ẏ = motion.uc[2](t)
    θ̇ = motion.uc[3](t)

    pos = [x, y]
    vel = [ẋ, ẏ]

    Rx = [cos(θ)  -sin(θ);  sin(θ)  cos(θ)]
    Ru = θ̇*[-sin(θ)  -cos(θ); cos(θ)  -sin(θ)]

    return TangentSE2(pos, vel, Rx, Ru)
end


"Type (3): the main driver functions that implement prescribed body motion"
function move_body!(body::RigidBody, t::Float64)
    move_body!(MotionType(body), body, t)
end

function move_body!(
    ::Type{MotionFunction},
    body::RigidBody,
    t::Float64
    )
    se2 = get_transformation(body.motion, t)
    for i=1:size(body.xb, 1)
        @views body.xb[i, :] .= se2.pos .+ se2.Rx*body.x0[i, :]
        @views body.ub[i, :] .= se2.vel .+ se2.Ru*body.x0[i, :]
    end
end
