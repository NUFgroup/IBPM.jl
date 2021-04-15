"""
From a method-of-lines perspective, this should contain everything related
    to the continuous-time problem... analogous to defining the RHS of an ODE.
    Then everything related to time-stepping goes in the IBProblem

At the moment, the distinction between this and IBProblem is not useful,
    but could potentially be useful for interfacing with DifferentialEquations
"""

"""
    IBMatrices(grid, bodies, Re)

Matrices that can be precomputed

Correspondence with Taira & Colonius (2007)
    Note that not all matrices defined here are explicitly constructed
C  - Basic curl operator for single-grid
        Call curl! function for multigrid to take into account boundary conditions
R  - Transforms velocity flux to circulation: gamma = R*q
G  - Discrete gradient operator
D  - Discrete divergence operator... D = -G'
E  - Maps fluxes to body motion, i.e. u_B = E*q
        Note that H = -E' is the regularization operator
A  - Implicit time-stepping operator for velocity flux
        A = I - (dt/2/h^2)*Lap
Q  - Q = [G E'] Averaging operator (used in the nonlinear term)

NOTE: Most of the functionality for this is in fluid-operators/lin.jl
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
    function IBMatrices(grid::T, bodies::Array{V, 1}, Re::Float64) where T <: Grid where V <: Body
        mats = new()
        mats.C = get_C(grid)  # Basic single-grid curl operator
        mats.Lap = mats.C'*mats.C/Re    # Laplacian
        mats.Λ = Lap_eigs(grid)

        # Plan DST
        # TODO: DO YOU STILL NEED THIS AFTER CREATING OPERATORS???
        mats.dst_plan = get_dst_plan(ones(Float64, grid.nx-1, grid.ny-1));
        mats.RCinv = get_lap_inv(grid, mats.Λ, mats.dst_plan);

        mats.Q = get_Q( grid );  # Averaging operator for nonlinear term
        mats.W = get_W( grid );  # WHAT IS THIS???

        mats.E = coupling_mat( grid, bodies )   # interface-coupling/ib_coupling.jl
        mats.RET = (mats.E*mats.C)'     # Precompute this mat-mat product

        return mats
    end
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
    function WorkingMemory(grid::Grid)
        mg = (grid isa UniformGrid) ? 1 : grid.mg  # Number of grid levels
        work = new()
        work.q1 = zeros(grid.nq, mg)
        work.q2 = zeros(grid.nq, mg)
        work.q3 = zeros(grid.nq, mg)
        work.q4 = zeros(grid.nq, mg)
        work.Γ1 = zeros(grid.nΓ, mg)
        work.Γ2 = zeros(grid.nΓ, mg)
        work.Γ3 = zeros(grid.nΓ, mg)
        work.bc = zeros(grid.nΓ)
        return work
    end
end


"""
SolnModel contains information about simulation parameters and stores
all static (non-time varying) matrices.
"""
abstract type SolnModel end

struct IBModel{T <: Grid, V <: Body} <: SolnModel
    grid::T
    bodies::Array{V, 1}         # Array of bodies
    Re::Float64                 # Reynolds number
    Uinf::Float64               # Free-stream velocity
    α::Float64                  # Angle of attack
    mats::IBMatrices            # Various precomputed sparse matrices
    XX::Union{Array{Float64, 2}, Nothing}  #  x-locations for computing rotational fluxes
    YY::Union{Array{Float64, 2}, Nothing}  #  y-locations for computing rotational fluxes
    function IBModel(grid::T,
                     bodies::Array{V, 1},
                     Re::Float64;
                     Uinf=1.0,
                     α=0.0,
                     xc=0.0,
                     yc=0.0) where {T <: Grid, V <: Body}
        mats = IBMatrices(grid, bodies, Re)

        # TODO: Put in different function??
        "Precompute grid locations for cross products (used for rotational flux)"
        if MotionType(bodies) == MovingGrid
            # Define x-xc, y-yc on all grids
            @assert length(bodies) == 1
            motion = bodies[1].motion
            mg = (grid isa UniformGrid) ? 1 : grid.mg  # Number of grid levels
            nx, ny, h = grid.nx, grid.ny, grid.h
            XX = zeros(nx*(ny-1), mg)  # Number of y-fluxes
            YY = zeros(ny*(nx-1), mg)  # Number of x-fluxes

            for lev=1:grid.mg
                hc = h*2^(lev-1);  # Coarse grid spacing

                ### y-coordinates for calculating x-fluxes
                y = @. ((1:ny)-0.5-ny/2)*hc + ny/2*h - grid.offy
                YY[:, lev] = (ones(nx-1)*(y.-yc)')[:]

                ### x-coordinates for calculating y-fluxes
                x = @. ((1:nx)-0.5-nx/2)*hc + nx/2*h - grid.offx
                XX[:, lev] = ((x.-xc)*ones(ny-1)')[:]
            end
        else
            XX, YY = nothing, nothing
        end
        return new{T, V}(grid, bodies, Re, Uinf, α, mats, XX, YY)
    end
end
