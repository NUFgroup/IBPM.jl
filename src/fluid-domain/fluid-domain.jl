"""
Types of flow grids

UniformGrid can probably be eliminated once the code is fairly stable

MOVE TO fluid-domain??
"""
abstract type Grid end

struct UniformGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    nu::Int
    mg::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    u_idx::Any
    v_idx::Any
    ω_idx::Any
    split_flux::Any
end

struct MultiGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    nu::Int
    mg::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    u_idx::Any
    v_idx::Any
    ω_idx::Any
    split_flux::Any
    LEFT::Int
    RIGHT::Int
    BOT::Int
    TOP::Int
end

function make_grid(nx::Int, ny::Int, offx::Float64, offy::Float64, len::Float64; mg=1::Int)
    nΓ  = (nx-1)*(ny-1)  # Number of circulation points
    nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
    nq = nu + nv;  # Total num of vel (flux) points
    h = len / nx;  # Grid spacing

    # Define indexing functions
    u(x_idx, y_idx) = @. x_idx + (nx+1)*(y_idx-1)'
    v(x_idx, y_idx) = @. x_idx + nx*(y_idx-1)' + nu
    ω(x_idx, y_idx) = @. x_idx + (nx-1)*(y_idx-1)'

    "Return a view to 2D arrays of fluxes"
    split_flux(q) = reshape(@view(q[1:nu]), nx+1, ny), reshape(@view(q[nu+1:end]), nx, ny+1)
    split_flux(q, lev) = reshape(@view(q[1:nu, lev]), nx+1, ny), reshape(@view(q[nu+1:end, lev]), nx, ny+1)

    if mg==1
        return UniformGrid(nx, ny, nΓ, nq, nu, mg, offx, offy, len, h,
            u, v, ω, split_flux)
    else
        # Predefine constant offsets for indexing boundary conditions
        left = 0;  right = ny+1
        bot = 2*(ny+1); top = 2*(ny+1) + nx+1
        return MultiGrid(nx, ny, nΓ, nq, nu, mg, offx, offy, len, h,
            u, v, ω, split_flux, left, right, bot, top)
    end

end
