
"""
function get_mats(grid::T, bodies::Array{V, 1}, Re::Float64) where T <: Grid where V <: Body
    return get_mats( MotionType(bodies), grid, bodies, Re)
end
"""

function get_mats(grid::T, bodies::Array{V, 1}, Re::Float64) where T <: Grid where V <: Body
    """
    Initialize all the matrices needed

    For Static motions, we can precompute E, E', and R*E'

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
    Q  - Q = [G E'] Averaging operator??
    """
    C = get_C(grid)
    Lap = C'*C/Re  # Laplacian

    Λ = Lap_eigs(grid)

    #  DOCUMENT THESE
    Q = get_Q( grid );
    W = get_W( grid );

    # Interpolation and regularization matrices
    E = coupling_mat( grid, bodies )   # ib_coupling.jl

    # Plan DST
    dst_plans = get_dst_plan(ones(Float64, grid.nx-1, grid.ny-1));

    lap_inv = get_lap_inv(grid, Λ, dst_plans);

    return IBMatrices(C, Lap, Λ, lap_inv, Q, W, E, (E*C)', dst_plans)
end


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


function update_coupling!(model::IBModel, t::Float64)
    bodies, grid = model.bodies, model.grid
    for j=1:length(bodies)
        move_body!(bodies[j].xb, bodies[j].ub, bodies[j].motion, t)
    end

    model.mats.E = coupling_mat( grid, bodies )
    model.mats.RET = (model.mats.E*model.mats.C)'
end

"""
For dispatching to Static motions

function static_fn(model::IBModel{<:Grid, <:Body{Static}})

For other Motions
function dynamic_fn(model::IBModel)
"""
