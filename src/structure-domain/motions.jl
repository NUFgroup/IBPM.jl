"""
Types of motion a body can undergo
"""
abstract type Motion end

abstract type InertialMotion <: Motion end

"Fixed bodies - no motion"
struct Static <: InertialMotion
end

"Rotating cylinder-specific type"
struct RotatingCyl <: InertialMotion
    Ω::Float64               # Angular velocity (constant)
end

"""
Body-fixed (non-inertial frame)... but really inertial equations on a moving grid

Express arbitrary motion as a combination of translation and a rotation
 => becomes a superposition of time-varying potential flow and
    solid-body rotation in base flux (no additional nonlinear terms needed)
See Hsieh-Chen Tsai thesis (2016) for derivation

Specify the linear and angular velocities relative to the lab frame

Note that the angular velocity θ is the negative of aerodynamic pitch.
This will override the definition of Uinf in the parent IBModel.
"""
mutable struct MovingGrid <: Motion
    U::Any            # Free-stream velocity U(t)
    θ::Any            # Angular position
    θ̇::Any            # Angular velocity
    xc::Float64       # x-center of rotation
    yc::Float64       # y-center of rotation
    MovingGrid(U, θ, θ̇; xc=0.0, yc=0.0) = new(U, θ, θ̇, xc, yc)
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

mutable struct MotionFunction <: InertialMotion
    xc::Any                  # Center position [x, y, θ] = xc(t)
    uc::Any                  # Center velocity [ẋ, ẏ, θ̇] = uc(t)
end
