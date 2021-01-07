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


"""
Motion type
"""
abstract type Motion end

# Fixed bodies - no motion
struct Static <: Motion end

"""
Body types
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T              # Motion function
    xb::Array{Float64, 2}  # (x, y) locations of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end


function get_body_info( bodies::Array{V, 1} where V <: Body )
    # determine the num of body points per body
    nf = [length(bodies[j].xb) for j=1:length(bodies)]
    nb = [size(bodies[j].xb, 1) for j=1:length(bodies)]
    return nb, nf
end

# Example body constructor
function make_cylinder(r, h, y0, motion)
    # build cylinder of radius r using flow grid spacing of h
    # (cylinder will have a spacing of 2h)

    circum = 2 * π * r; #  Circumference of the circle

    # Get # of points such that ds = 2h
    n = Int( floor( circum / h / 2 ) );

    int =  2*π/n ;
    spt = 0:int:(n-1)*int;
    xhat = r.*cos.(spt);
    yhat = r.*sin.(spt);

    xb = [xhat  yhat.+y0];

    # sanity check: make sure ds is equal to 2 * h
    ds = sqrt( (xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2 ) ;

    return RigidBody(motion, xb, fill(ds, n))
end




function make_naca(x0, N, spec, motion)

    # Define x-locations
    dθ = π/(N-1)

    x = 0.5.*(1 .+ cos.(dθ.*(0:N-1)))
    _, xU, xL, yU, yL = naca(x, spec);

    # Edge points of "panels"
    xe = zeros(2*N-1, 2)
    for i=1:N
        xe[i, :] = [x[i]-x0, yU[i]]
    end
    for i=1:N-1
        xe[i+N, :] = [x[N-i]-x0, yL[N-i]]
    end

    # What we'll actually keep is the center points and edge lengths
    #   (as in vortex panel methods)
    xb = (xe[2:end, :] .+ xe[1:end-1, :] ) ./ 2.0
    ds = sqrt.( sum( (xe[2:end, :] .- xe[1:end-1, :]).^2, dims=2) )

    return RigidBody(motion, xb, ds[:, 1])
end
