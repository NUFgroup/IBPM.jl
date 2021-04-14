"""
Body types

See structure-domain/sample-bodies.jl for examples of constructing these

MOVE TO structure-domain???
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T               # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Reference locations (for moving bodies)
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds::Array{Float64, 1}   # line segment lengths on body
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
