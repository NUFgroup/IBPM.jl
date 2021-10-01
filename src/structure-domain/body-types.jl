"""
Body types

See structure-domain/sample-bodies.jl for examples of constructing these
"""
abstract type Body{T <: Motion} end

struct RigidBody{T <: Motion} <: Body{T}
    motion::T              # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Reference locations (for moving bodies)
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end

struct EB_BCinfo
    BCvec :: Vector{Any}
    BCnode::Float64
end

struct DeformingBody{T<: Motion} <: Body{T}
    motion::T              # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Original locations (for moving bodies)
    xref::Array{Float64, 2}  # Reference locations (about which displacements
                             # are determined)
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds0::Array{Float64, 1}   # line segment lengths on body (in undeformed configuration)
    ds::Array{Float64, 1}   # line segment lengths on body (in deformed configuration)
    kb::Vector{Float64} #array of structural bending stiffness values
    ke::Vector{Float64} #array of structural extensional stiffness values
    m::Vector{Float64} #array of structural mass values
    BCinfo::Array{EB_BCinfo,1}
end
