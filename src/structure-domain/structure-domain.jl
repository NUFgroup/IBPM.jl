
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


function get_transformation(
        motion::MotionFunction,
        t::Float64   # Should be time of the *next* step (i.e. t[k+1])
        )
    x, y, θ = motion.xc(t)
    ẋ, ẏ, θ̇ = motion.uc(t)

    pos = [x, y]
    vel = [ẋ, ẏ]

    Rx = [cos(θ)  -sin(θ);  sin(θ)  cos(θ)]
    Ru = θ̇*[-sin(θ)  -cos(θ); cos(θ)  -sin(θ)]

    return TangentSE2(pos, vel, Rx, Ru)
end


function move_body!(body::RigidBody, t::Float64)
    move_body!(MotionType(body), body, t)
end

function move_body!(
    ::Type{RotatingCyl},
    body::RigidBody,
    t::Float64
    )
    """
    Compute linear speed of points on rotating cylinder
    θ = atan(y/x)
    ẋ = -R*Ω*sin(θ)
    ẏ =  R*Ω*cos(θ)
    Note: this is a specialized case for development only
    """
    xb = body.xb
    ub = body.ub
    motion = body.motion
    R = sqrt(xb[1, 1]^2 + xb[1, 2]^2)
    θ = atan.(xb[:, 2], xb[:, 1])
    ub[:, 1] .= -R*motion.Ω*sin.(θ)
    ub[:, 2] .=  R*motion.Ω*cos.(θ)
    return ub
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


# TODO: DO WE NEED THIS?? SEEMS LIKE MAYBE NOT...
function move_body!(
    ::Type{BodyFixed},
    body::RigidBody,
    t::Float64
    )
    body.ub .*= 0.0
end
