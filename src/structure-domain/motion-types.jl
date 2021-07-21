"""
Types of motion a body can undergo
"""
abstract type Motion end

abstract type InertialMotion <: Motion end #for cases where grid is not moving

"Fixed body - no motion"
struct Static <: InertialMotion
end

"Body moves with prescribed kinematics in the domain"
#Case 1: user prescribes body motion via MotionFunction
mutable struct MotionFunction <: InertialMotion
    xc::Any                  # Center position [x, y, θ] = xc(t)
    uc::Any                  # Center velocity [ẋ, ẏ, θ̇] = uc(t)
end

#Deprecated??
# "Rotating cylinder-specific type"
# struct RotatingCyl <: InertialMotion
#     Ω::Float64               # Angular velocity (constant)
# end

#Case 2: user prescribes body motion by moving the grid and adding the necessary
#   inertial terms.
# **Preferred where possible, leads to MUCH faster computations because
# expensive operators can be pre-computed when the body is treated as stationary
# in the inertial frame.
#
# The approach here expresses arbitrary motion as a combination of translation
# and a rotation
#  => becomes a superposition of time-varying potential flow and
#     solid-body rotation in base flux (no additional nonlinear terms needed)
#
#     See Hsieh-Chen Tsai thesis (2016) for derivation
#
# Specify the linear and angular velocities relative to the lab frame
#
# Note that the angular velocity θ is the negative of aerodynamic pitch.

mutable struct MovingGrid <: Motion
    U::Any            # Free-stream x-velocity U(t)
    V::Any            # Free-stream y-velocity V(t)
    θ::Any            # Angular position (t)
    θ̇::Any            # Angular velocity (t)
    xc::Float64       # x-center of rotation
    yc::Float64       # y-center of rotation
    MovingGrid(U, V, θ, θ̇; xc=0.0, yc=0.0) = new(U, V, θ, θ̇, xc, yc)
end
