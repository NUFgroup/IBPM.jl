
"""
Motion type
"""
abstract type Motion end

# Fixed bodies - no motion
struct Static <: Motion end

# Rotating cylinder-specific type
struct RotatingCyl <: Motion
    θ̇::Float64        # Angular velocity
end


function move_body!(
    xb::Array{Float64, 2},
    ub::Array{Float64, 2},
    motion::RotatingCyl,
    t::Float64
    )
    """
    Compute linear speed of points on rotating cylinder
    θ = atan(y/x)
    ẋ = -R*θ̇*sin(θ)
    ẏ =  R*θ̇*cos(θ)

    Note: this is a specialized case for development only
    """
    R = sqrt(xb[1, 1]^2 + xb[1, 2]^2)
    θ = atan.(xb[:, 2], xb[:, 1])
    ub[:, 1] .= -R*motion.θ̇*sin.(θ)
    ub[:, 2] .=  R*motion.θ̇*cos.(θ)
    return ub
end


"""
Struct to hold maps for position and velocity
    Modeled after Rowley's C++ code

    position transformation:
        x = pos .+ Rx*x
    velocity transformation:
        v = vel .+ Rv*v
"""
struct TangentSE2
    pos::Array{Float64, 1}   # [x, y] position of center
    vel::Array{Float64, 1}   # [ẋ, ẏ] velocity of center
    Rx::Any  # maps positions [x, y]
    Ru::Any  # maps velocities [ẋ, ẏ]
end


function move_body!(
    xb::Array{Float64, 2},
    ub::Array{Float64, 2},
    motion::Motion,
    t::Float64
    )
    se2 = get_transformation(motion, t)
    #println(ub[1, 1])
    #println(se2.Ru[1, :] * ub[1, :])
    for i=1:size(xb, 1)
        @views xb[i, :] .= se2.pos .+ se2.Rx*xb[i, :]
        @views ub[i, :] .= se2.vel .+ se2.Ru*xb[i, :]
    end
    #println(se2.Ru)
    return ub
    #xb .= se2.pos .+ se2.Rx*xb
    #ub .= se2.vel .+ se2.Ru*ub
end

mutable struct MotionFunction <: Motion
    xc::Any                  # Center position [x, y, θ] = xc(t)
    uc::Any                  # Center velocity [ẋ, ẏ, θ̇] = uc(t)
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
