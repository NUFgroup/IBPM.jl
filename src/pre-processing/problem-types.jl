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

abstract type InertialMotion <: Motion end

# Fixed bodies - no motion
struct Static <: InertialMotion
end

# Rotating cylinder-specific type
struct RotatingCyl <: InertialMotion
    Ω::Float64               # Angular velocity (constant)
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

"""
Body-fixed (non-inertial) reference frame

Express arbitrary motion as a combination of translation and a rotation
 => becomes a superposition of time-varying potential flow and
    solid-body rotation in base flux (no additional nonlinear terms needed)
See Hsieh-Chen Tsai thesis (2016) for derivation

Specify the linear and angular velocities relative to the lab frame

Note that the angular velocity θ is the negative of aerodynamic pitch.
This will override the definition of Uinf in the parent IBModel
"""
mutable struct BodyFixed <: Motion
    U::Any            # Free-stream velocity U(t)
    θ::Any            # Angular position
    θ̇::Any            # Angular velocity
end

"""
Body types
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T               # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Reference locations (for moving bodies)
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end

"""
Matrices that can be precomputed
"""

mutable struct IBMatrices
    C::SparseArrays.SparseMatrixCSC{Float64,Int64}
    Lap::SparseArrays.SparseMatrixCSC{Float64,Int64}
    Λ::Array{Float64,2}
    RCinv::LinearMap
    Q::SparseArrays.SparseMatrixCSC{Float64,Int64}
    W::SparseArrays.SparseMatrixCSC{Float64,Int64}
    E::AbstractArray
    RET::AbstractArray
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
    Uinf::Float64               # Free-stream velocity
    α::Float64                  # Angle of attack
    mats::IBMatrices            # Various precomputed sparse matrices
end

"""
IBProblem has the info of IBModel as well as the problem structure (e.g., the
explicit time stepping scheme and information about the implicit treatment via
the A and B matrices and their inverses)
"""
abstract type AbstractIBProblem end

mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    work::WorkingMemory
    A
    Ainv
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
