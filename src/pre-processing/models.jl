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
C' - Transforms velocity flux to circulation: Γ = C'*q
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
    C::LinearMap
    Λ::Array{Float64, 2}
    Δinv::LinearMap
    E::LinearMap
    #Q::SparseArrays.SparseMatrixCSC{Float64,Int64}
    #W::SparseArrays.SparseMatrixCSC{Float64,Int64}
    dst_plan::Any
    function IBMatrices(grid::T, bodies::Array{V, 1}) where T <: Grid where V <: Body
        mats = new()
        mats.C = LinearMap( (q, ψ) -> ibpm.curl!(q, ψ, grid),  # Forward
                            (Γ, q) -> ibpm.rot!(Γ, q, grid),   # Transpose
                            grid.nq, grid.nΓ; ismutating=true)
        #mats.Lap = mats.C'*mats.C/Re    # Laplacian
        mats.Λ = lap_eigs(grid)

        # Plan DST for inverse Laplacian
        mats.dst_plan = get_dst_plan(ones(Float64, grid.nx-1, grid.ny-1));
        mats.Δinv = get_lap_inv(grid, mats.Λ, mats.dst_plan)

        # Interpolation/regularization matrix
        mats.E = ibpm.setup_reg(grid, bodies)   # interface-coupling/interface-oupling.jl

        # Averaging operators for nonlinear term
        #mats.Q = get_Q( grid );
        #mats.W = get_W( grid );
        return mats
    end
end

"""
Different time stepping schemes
"""
abstract type ExplicitScheme end

# TODO: Make constructor to generate β automatically
struct AdamsBashforth <: ExplicitScheme
    dt::Float64
    β::Array{Float64, 1}
end

"""
Pre-allocate memory to certain vectors that can be re-used throughout the
computation process
"""
mutable struct WorkingMemory
    q1::Array{Float64, 2}
    q2::Array{Float64, 2}
    q3::Array{Float64, 2}
    q4::Array{Float64, 2}
    q5::Array{Float64, 2}
    Γ1::Array{Float64, 2}
    Γ2::Array{Float64, 2}
    Γ3::Array{Float64, 2}
    Γbc::Array{Float64, 1}    # Poisson boundary conditions for multigrid
    rhsbc::Array{Float64, 1}  # Time-stepping boundary conditions for multigrid
    function WorkingMemory(grid::Grid)
        work = new()
        work.q1 = zeros(grid.nq, grid.mg)
        work.q2 = zeros(grid.nq, grid.mg)
        work.q3 = zeros(grid.nq, grid.mg)
        work.q4 = zeros(grid.nq, grid.mg)
        work.q5 = zeros(grid.nq, grid.mg)
        work.Γ1 = zeros(grid.nΓ, grid.mg)
        work.Γ2 = zeros(grid.nΓ, grid.mg)
        work.Γ3 = zeros(grid.nΓ, grid.mg)
        work.Γbc = zeros(2*(grid.nx+1)+2*(grid.ny+1))
        work.rhsbc = zeros(grid.nΓ)
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
    work::WorkingMemory
    XX::Union{Array{Float64, 2}, Nothing}  #  x-locations for computing rotational fluxes
    YY::Union{Array{Float64, 2}, Nothing}  #  y-locations for computing rotational fluxes
    function IBModel(grid::T,
                     bodies::Array{V, 1},
                     Re::Number;
                     Uinf=1.0,
                     α=0.0,
                     xc=0.0,
                     yc=0.0) where {T <: Grid, V <: Body}
        mats = IBMatrices(grid, bodies)
        work = WorkingMemory(grid)

        # TODO: Put in different function??
        "Precompute grid locations for cross products (used for rotational flux)"
        if MotionType(bodies) == MovingGrid
            # Define x-xc, y-yc on all grids
            @assert length(bodies) == 1
            motion = bodies[1].motion
            nx, ny, h = grid.nx, grid.ny, grid.h
            XX = zeros(nx*(ny+1), grid.mg)  # Number of y-fluxes
            YY = zeros(ny*(nx+1), grid.mg)  # Number of x-fluxes

            for lev=1:grid.mg
                hc = h*2^(lev-1);  # Coarse grid spacing

                ### y-coordinates for calculating x-fluxes
                y = @. ((1:ny)-0.5-ny/2)*hc + ny/2*h - grid.offy
                YY[:, lev] = (ones(nx+1)*(y.-yc)')[:]

                ### x-coordinates for calculating y-fluxes
                x = @. ((1:nx)-0.5-nx/2)*hc + nx/2*h - grid.offx
                XX[:, lev] = ((x.-xc)*ones(ny+1)')[:]
            end
        else
            XX, YY = nothing, nothing
        end
        return new{T, V}(grid, bodies, Float64(Re), Uinf, α, mats, work, XX, YY)
    end
end
