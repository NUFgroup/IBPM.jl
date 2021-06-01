"""
Functions used in the pre-processing stage, from initializing the state
variables and problem structure to building and storing all static (non-time
varying) matrices
"""


"""
initialize the state variables (stores everything needed for time stepping)
"""
function init_state(prob::IBProblem)
        grid = prob.model.grid
        mg = (grid isa UniformGrid) ? 1 : grid.mg  # Number of grid levels

        nb, nf = get_body_info(prob.model.bodies)
        return IBState{typeof(grid)}(
        zeros(grid.nq, mg),   # Flux
        zeros(grid.nq, mg),   # Background flux
        zeros(grid.nΓ, mg),   # Circulation
        zeros(grid.nΓ, mg),   # Streamfunction
        [zeros(grid.nΓ, mg) for i=1:length(prob.scheme.β)],   # Nonlinear
        [zeros(nf[i]) for i=1:length(nf)],   # Surface stresses on each body
        zeros(sum(nf)),
        [zeros(length(prob.model.bodies))';],
        zeros(length(prob.model.bodies)),
        0.0,
        0.0
        )
end

"""
initialize the simulation parameters and static matrices
"""
function init_model(grid::T where T <: Grid,
                    bodies::Array{V, 1} where V <: Body,
                    Re::Float64)
    return IBModel(grid, bodies, Re, get_mats(grid, bodies, Re))

end


"""
initialize the problem structure (matrices used, bodies and simulation
parameters, time steppping scheme, ...)

Note: the scheme actually speaks to the terms that are explicitly treated. This
is a projection method the directly enforces the no-slip condition, so some terms
are implicitly treated. This information is not contained in scheme, but in the
A, Ainv, B, and Binv matrices
"""
function init_prob(grid::T where T <: Grid,
                   bodies::Array{V, 1} where V <: Body,
                   Re::Float64,
                   dt::Float64,
                   freestream::NamedTuple)
    mats = get_mats(grid, bodies, Re)
    model = IBModel(grid, bodies, Re, mats, freestream)
    scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
    A, Ainv, Binv = get_AB(model, dt)
    work = init_memory(grid)
    return IBProblem(model, scheme, work, A, Ainv, Binv)
end


"""
For optimal performance, pre-allocate a few vectors that can be re-used
throughout the code
"""

function init_memory(grid::UniformGrid)
    return WorkingMemory(
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ)
    )
end

function init_memory(grid::MultiGrid)
    return WorkingMemory(
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ)
    )
end

"""
Build the matrices that can be pre-computed because they are time-constant
"""

function get_mats(grid::T, bodies::Array{V, 1}, Re::Float64) where T <: Grid where V <: Body
    """
    Initialize all the matrices needed

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
    """
    C = get_C(grid)
    Lap = C'*C/Re  # Laplacian

    Λ = Lap_eigs(grid)

    Q = get_Q( grid );
    W = get_W( grid );

    E = coupling_mat( grid, bodies )   # ib_coupling.jl

    # Plan DST
    #    DO YOU STILL NEED THESE AFTER CREATING OPERATORS???
    dst_plans = get_dst_plan(ones(Float64, grid.nx-1, grid.ny-1));

    lap_inv = get_lap_inv(grid, Λ, dst_plans);
    return IBMatrices(C, Lap, Λ, lap_inv, Q, W, E, (E*C)', dst_plans)
end

"""
Compute the matrix (as a LinearMap) that represents the modified Poisson
operator (I + dt/2 * Beta * RC) arising from the implicit treatment of the
Laplacian. A system involving this matrix is solved to compute a trial
circulation that doesn't satisfy the BCs, and then again to use the surface
stresses to update the trial circulation so that it satisfies the BCs
"""
function get_A(model::IBModel{UniformGrid, <:Body}, dt::Float64)
    A = I - (dt/2 / (model.grid.h^2))*model.mats.Lap
    Ainv = get_Ainv(model, dt)
    return A, Ainv
end

function get_A(model::IBModel{MultiGrid, <:Body}, dt::Float64)
    hc = [model.grid.h * 2^(lev-1) for lev=1:model.grid.mg]
    A = [I - (dt/2 / (hc[lev]^2)) *model.mats.Lap for lev=1:model.grid.mg]
    Ainv = [get_Ainv(model, dt; lev=lev) for lev=1:model.grid.mg]
    return A, Ainv
end

"""
Pre-compute and store the matrix, B, that is used to solve for the surface
stresses that lead to a velocity that satisfies the no-slip BCs
"""
function get_B(model::IBModel{UniformGrid, RigidBody{Static}}, dt::Float64, Ainv::LinearMap)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # need to build and store surface stress matrix and its inverse if at first time step
    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot, 1 );         # Unit vector

    # TODO: Alternative... could create a dummy state to operate on here
    b = zeros( nftot, 1 );         # Working array
    Γ = zeros(model.grid.nΓ, 1)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, 1)    # Working array for streamfunction
    q = zeros(model.grid.nq, 1)    # Working array for velocity flux

    for j = 1 : nftot
        if j>1
            e[j-1] = 0.0
        end
        e[j] = 1.0;

        b_times!( b, e, Ainv, model, Γ, ψ, q );
        B[:, j] = b
    end
    Binv = inv(B)
    return B, Binv
end


function get_B(model::IBModel{MultiGrid, RigidBody{Static}}, dt::Float64, Ainv)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # need to build and store surface stress matrix and its inverse if at first time step
    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot, 1 );         # Unit vector

    # TODO: Alternative... could create a dummy state to operate on here
    b = zeros( nftot, 1 );         # Working array
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    for j = 1 : nftot
        if j>1
            e[j-1] = 0.0
        end
        e[j] = 1.0;

        # Only fine-grid Ainv is used here
        b_times!( b, e, Ainv[1], model, Γ, ψ, q );
        B[:, j] = b
    end
    Binv = inv(B)
    return B, Binv
end


"""
Create A and B
"""
function get_AB(model::IBModel, dt::Float64)
    A, Ainv = get_A(model, dt)
    B, Binv = get_B(model, dt, Ainv)
    return A, Ainv, B, Binv
end
