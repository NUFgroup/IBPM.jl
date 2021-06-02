"""
Types of flow grids

UniformGrid can probably be eliminated once the code is fairly stable

MOVE TO fluid-domain??
"""
abstract type Grid end

struct MultiGrid <: Grid
    nx::Int
    ny::Int
    nΓ::Int
    nq::Int
    mg::Int
    offx::Float64
    offy::Float64
    len::Float64
    h::Float64
    split_flux::Any
    LEFT::Int
    RIGHT::Int
    BOT::Int
    TOP::Int
end

# function make_grid(nx::Int, ny::Int, offx::Float64, offy::Float64, len::Float64; mg=1::Int)
#     nΓ  = (nx-1)*(ny-1)  # Number of circulation points
#     nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
#     nq = nu + nv;  # Total num of vel (flux) points
#     h = len / nx;  # Grid spacing
#
    # "Return views to 2D arrays of fluxes"
    # split_flux(q; lev=1) = reshape(@view(q[1:nu, lev]), nx+1, ny),
    #                        reshape(@view(q[nu+1:end, lev]), nx, ny+1)
    #
    # # Predefine constant offsets for indexing boundary conditions
    # left = 0;  right = ny+1
    # bot = 2*(ny+1); top = 2*(ny+1) + nx+1
    # return MultiGrid(nx, ny, nΓ, nq, mg, offx, offy, len, h,
    #     split_flux, left, right, bot, top)
# end

function make_grid(h::Float64, boundary::NTuple{4,Float64}; mg=1::Int)

    #back out variables that the software needs from user defined vars
    offx = -boundary[1]
    offy = -boundary[3]
    len = boundary[2]-boundary[1]
    ylen = boundary[4]-boundary[3]

    nx = Int64(round(len/h))
    ny = Int64(round(ylen*nx/len))

    nΓ  = (nx-1)*(ny-1)  # Number of circulation points

    nu = ny * (nx+1); nv = nx * (ny+1);  # num of (flux) points
    # Total num of vel (flux) points
    nq = nu + nv;

    h = len / nx;  # Grid spacing

    "Return views to 2D arrays of fluxes"
    split_flux(q; lev=1) = reshape(@view(q[1:nu, lev]), nx+1, ny),
                           reshape(@view(q[nu+1:end, lev]), nx, ny+1)

    # Predefine constant offsets for indexing boundary conditions
    left = 0;  right = ny+1
    bot = 2*(ny+1); top = 2*(ny+1) + nx+1

    return MultiGrid(nx, ny, nΓ, nq, mg, offx, offy, len, h,
        split_flux, left, right, bot, top)
end
