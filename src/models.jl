

function get_mats(grid::T, bodies::Array{V, 1}, Re::Float64) where T <: Grid where V <: Body
    """
    Initialize all the matrices needed

    Correspondence with Taira & Colonius (2007)
        Note that not all matrices defined here are explicitly constructed
    C  - Basic curl operator for single-grid
            Call curl! function for multigrid to take into account boundary conditions
    R  - Transforms velocity to velocity flux: q = R*u
    G  - Discrete gradient operator
    D  - Discrete divergence operator... D = -G'
    E  - Maps fluxes to body motion, i.e. u_B = E*q
            Note that H = -E' is the regularization operator
    A  - Implicit time-stepping operator for velocity flux
            A = I - (dt/2/h^2)*Lap
    Q  - Q = [G E'] Averaging operator
    """
    C = get_C(grid)
    R = C'        # Transforms velocity to flux, i.e. q = R*u
    Lap = R*C/Re  # Laplacian

    Λ = Lap_eigs(grid)

    Q = get_Q( grid );
    W = get_W( grid );

    ET = reg_mats( grid, bodies )   # ib_coupling.jl
    E = ET';

    # Plan DST
    #    DO YOU STILL NEED THESE AFTER CREATING OPERATORS???
    dst_plans = get_dst_plan(ones(Float64, grid.nx-1, grid.ny-1));

    RCinv = get_RCinv(grid, Λ, dst_plans);

    return IBMatrices(C, R, Lap, Λ, RCinv, Q, W, E, ET, R*ET, dst_plans)
end


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

abstract type NavierStokesModel end

# TODO Should probably be called something else...
struct IBModel{T <: Grid, V <: Body} <: NavierStokesModel
    grid::T
    bodies::Array{V, 1}         # Array of bodies
    Re::Float64                 # Reynolds number
    mats::IBMatrices            # Various precomputed sparse matrices
end

function init_model(grid::T where T <: Grid,
                    bodies::Array{V, 1} where V <: Body,
                    Re::Float64)
    return IBModel(grid, bodies, Re, get_mats(grid, bodies, Re))

end



"""
For dispatching to Static motions

function static_fn(model::IBModel{<:Grid, <:Body{Static}})

For other Motions
function dynamic_fn(model::IBModel)
"""
