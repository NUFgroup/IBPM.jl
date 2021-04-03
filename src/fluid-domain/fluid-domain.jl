
function make_grid(h::Float64, boundary::NTuple{4,Float64}; mg=1::Int)

    #back out variables that the software needs from user defined vars
    offx = -boundary[1]
    offy = -boundary[3]
    len = boundary[2]-boundary[1]
    ylen = boundary[4]-boundary[3]

    nx = Int64(round(len/h))
    ny = Int64(round(ylen*nx/len))

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
