

"""
Body types
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T               # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    ub::Array{Float64, 2}   # (ẋ, ẏ) velocities of body points
    ds::Array{Float64, 1}   # line segment lengths on body
end


function get_ub(bodies::Array{V, 1} where V<:Body)
    nb, nf = get_body_info(bodies)
    ub = zeros(Float64, sum(nb), 2)

    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j=1:length(bodies)
        #move_body!(bodies[j].xb, bodies[j].ub, bodies[j].motion, t)  # Not time-dependent
        ub[nbod_tally .+ (1:nb[j]), :] .= bodies[j].ub

        # update body index
        nbod_tally += nb[j];
    end
    return [ub[:, 1]; ub[:, 2]]  # Stack [ub; vb]
end


function get_body_info( bodies::Array{V, 1} where V <: Body )
    # determine the num of body points per body
    nf = [length(bodies[j].xb) for j=1:length(bodies)]
    nb = [size(bodies[j].xb, 1) for j=1:length(bodies)]
    return nb, nf
end

# Example body constructor
function make_cylinder(r, h, y0, motion; n=0)
    # build cylinder of radius r using flow grid spacing of h
    # (cylinder will have a spacing of 2h)

    circum = 2 * π * r; #  Circumference of the circle

    # Get # of points such that ds = 2h
    if (n==0)
        n = Int( floor( circum / h / 2 ) );
    end

    int =  2*π/n ;
    spt = 0:int:(n-1)*int;
    xhat = r.*cos.(spt);
    yhat = r.*sin.(spt);

    xb = [xhat  yhat.+y0];

    # sanity check: make sure ds is equal to 2 * h
    ds = sqrt( (xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2 ) ;

    return RigidBody(motion, xb, 0.0*xb, fill(ds, n))
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

    return RigidBody(motion, xb, 0.0*xb, ds[:, 1])
end


function MotionType( bodies::Array{V, 1} where V<:Body )
    motions = [typeof(bodies[i].motion) for i=1:length(bodies)]
    if all(motions .== Static)
        return Static
    elseif all(motions .== RotatingCyl)
        return RotatingCyl
    else
        return Motion
    end
end
