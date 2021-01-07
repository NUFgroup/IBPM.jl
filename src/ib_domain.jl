"""
Grid types
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


function make_grid(nx::Int, ny::Int, offx::Float64, offy::Float64, len::Float64; mg=1::Int)
    nΓ  = (nx-1)*(ny-1)  # Number of circulation points

    # num of (flux) points
    nu = ny * (nx-1); nv = nx * (ny-1);
    # Total num of vel (flux) points
    nq = nu + nv;

    h = len / nx;  # Grid spacing


    if mg==1
        return UniformGrid(nx, ny, nΓ, nq, offx, offy, len, h)
    else
        return MultiGrid(nx, ny, nΓ, nq, offx, offy, len, h, mg,
                         zeros(nΓ, 4)  # Streamfunction boundary conditions
                         )
    end
end
