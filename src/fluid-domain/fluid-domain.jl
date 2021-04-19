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
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    u_idx::Any
    v_idx::Any
    ω_idx::Any
end

struct MultiGrid <: Grid
end

function make_grid(nx::Int, ny::Int, offx::Float64, offy::Float64, len::Float64; mg=1::Int)
    nΓ  = (nx-1)*(ny-1)  # Number of circulation points
    nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
    nq = nu + nv;  # Total num of vel (flux) points
    h = len / nx;  # Grid spacing

    if mg==1
        # Define indexing functions
        # TODO: move to arrays?? This allocates memory
        u(x_idx, y_idx) = @. x_idx + (nx+1)*(y_idx-1)'
        v(x_idx, y_idx) = @. x_idx + nx*(y_idx-1)' + nu
        ω(x_idx, y_idx) = @. x_idx + (nx-1)*(y_idx-1)'
        #u = (1:nx+1) .+ (nx+1)*(0:ny-1)'
        #v = (1:nx) .+ nx*(0:ny)' .+ nu
        #ω = (1:nx-1) .+ (nx-1)*(0:ny-2)'
        return UniformGrid(nx, ny, nΓ, nq, offx, offy, len, h, u, v, ω)
    end
end
