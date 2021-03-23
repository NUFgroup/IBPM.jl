
function get_body_info( bodies::Array{V, 1} where V <: Body )
    # determine the num of body points per body
    nf = [length(bodies[j].xb) for j=1:length(bodies)]
    nb = [size(bodies[j].xb, 1) for j=1:length(bodies)]
    return nb, nf
end


function MotionType(prob::IBProblem)
    bodies = prob.model.bodies
    motions = [typeof(bodies[i].motion) for i=1:length(bodies)]
    if all(motions .== Static)
        return Static
    else
        return Motion
    end
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



#Why is this here?
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
