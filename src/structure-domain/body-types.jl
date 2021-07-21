"""
Body types

See structure-domain/sample-bodies.jl for examples of constructing these

MOVE TO structure-domain???
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T              # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Reference locations (for moving bodies)
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end
