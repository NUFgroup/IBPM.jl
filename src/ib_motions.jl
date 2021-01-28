
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


mutable struct MotionFunction <: Motion
    xc::Any                  # Center position [x, y, θ] = xc(t)
    uc::Any                  # Center velocity [ẋ, ẏ, θ̇] = uc(t)
end
