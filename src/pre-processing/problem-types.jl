"""
Types and structs used to define cross-disciplinary variables

Maybe some opportunity for restructuring...
"""

"""
Types of flow grids
"""
abstract type Grid end

struct UniformGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
end

struct MultiGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    mg::Int
    stbc::Array{Float64, 2}
end


"""
Types of motion a body can undergo
"""
abstract type Motion end

# Fixed bodies - no motion
struct Static <: Motion end

"""
Body types
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T              # Motion function
    xb::Array{Float64, 2}  # (x, y) locations of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end

"""
Matrices that can be precomputed
"""
struct IBMatrices
    C::SparseArrays.SparseMatrixCSC{Float64,Int64}
    R::SparseArrays.SparseMatrixCSC{Float64,Int64}
    Lap::SparseArrays.SparseMatrixCSC{Float64,Int64}
    Λ::Array{Float64,2}
    RCinv::LinearMap
    Q::SparseArrays.SparseMatrixCSC{Float64,Int64}
    W::SparseArrays.SparseMatrixCSC{Float64,Int64}
    E::SparseArrays.SparseMatrixCSC{Float64,Int64}
    ET::SparseArrays.SparseMatrixCSC{Float64,Int64}
    RET::SparseArrays.SparseMatrixCSC{Float64,Int64}
    dst_plan::Tuple{Any, Array{Float64, 2}}
end

"""
Different time stepping schemes
"""
abstract type ExplicitScheme end


struct AdamsBashforth <: ExplicitScheme
    dt::Float64
    β::Array{Float64, 1}
end

"""
Pre-allocate memory to certain vectors that can be re-used throughout the
computation process
"""
mutable struct WorkingMemory
    q1::AbstractArray
    q2::AbstractArray
    q3::AbstractArray
    q4::AbstractArray
    Γ1::AbstractArray
    Γ2::AbstractArray
    Γ3::AbstractArray
    bc::AbstractArray
end


"""
SolnModel contains information about simulation parameters and stores
all static (non-time varying) matrices
"""
abstract type SolnModel end

struct IBModel{T <: Grid, V <: Body} <: SolnModel
    grid::T
    bodies::Array{V, 1}         # Array of bodies
    Re::Float64                 # Reynolds number
    mats::IBMatrices            # Various precomputed sparse matrices
end

"""
IBProblem has the info of IBModel as well as the problem structure (e.g., the
explicit time stepping scheme and information about the implicit treatment via
the A and B matrices and their inverses)
"""
struct IBProblem
    model::IBModel
    scheme::ExplicitScheme
    work::WorkingMemory
    A
    Ainv
    B
    Binv
end

"""
state variables (stores everything needed for time stepping)
"""
abstract type State end

mutable struct IBState{T<:Grid} <: State
        q::Array{Float64, 2}
        q0::Array{Float64, 2}
        Γ::Array{Float64, 2}     # Circulation
        ψ::Array{Float64, 2}     # Streamfunction
        nonlin::Array{Array{Float64, 2}, 1}  # Memory of nonlinear terms
        fb::Array{Array{Float64, 1}, 1}          # Surface stresses
        F̃b::Array{Float64, 1}                    # Body forces * dt
        CD::Array{Float64, 1}    # Drag coefficient
        CL::Array{Float64, 1}    # Lift coefficient
        cfl::Float64
        slip::Float64
end
