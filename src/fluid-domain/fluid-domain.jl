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
    u::Any
    v::Any
    ω::Any
end

function make_grid(nx::Int, ny::Int, offx::Float64, offy::Float64, len::Float64; mg=1::Int)
    nΓ  = (nx-1)*(ny-1)  # Number of circulation points
    nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
    nq = nu + nv;  # Total num of vel (flux) points
    h = len / nx;  # Grid spacing

    if mg==1
        u(i, j) = (nx+1)*(j-1) + i
        v(i, j) = nu + nx*(j-1) + i
        ω(i, j) = (ny-1)*(j-1) + i
        return UniformGrid(nx, ny, nΓ, nq, offx, offy, len, h, u, v, ω)
    end
end
