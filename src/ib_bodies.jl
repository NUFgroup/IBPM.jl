

"""
Body types
"""
abstract type Body{T <: Motion} end

struct RigidBody{T} <: Body{T}
    motion::T               # Motion function
    xb::Array{Float64, 2}   # (x, y) locations of body points
    x0::Array{Float64, 2}   # Reference locations (for moving bodies)
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
function make_cylinder(r, h, x0, y0, motion; n=0)
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

    xb = [xhat.+x0  yhat.+y0];

    # sanity check: make sure ds is equal to 2 * h
    ds = sqrt( (xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2 ) ;

    return RigidBody(motion, xb, copy(xb), 0.0*xb, fill(ds, n))
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

    return RigidBody(motion, xb, copy(xb), 0.0*xb, ds[:, 1])
end


function MotionType( body::V where V<:Body )
    return typeof(body.motion)
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


function get_transformation(
        motion::MotionFunction,
        t::Float64   # Should be time of the *next* step (i.e. t[k+1])
        )
    x, y, θ = motion.xc(t)
    ẋ, ẏ, θ̇ = motion.uc(t)

    pos = [x, y]
    vel = [ẋ, ẏ]

    Rx = [cos(θ)  -sin(θ);  sin(θ)  cos(θ)]
    Ru = θ̇*[-sin(θ)  -cos(θ); cos(θ)  -sin(θ)]

    return TangentSE2(pos, vel, Rx, Ru)
end


function move_body!(body::RigidBody, t::Float64)
    move_body!(MotionType(body), body, t)
end


function move_body!(
    ::Type{RotatingCyl},
    body::RigidBody,
    t::Float64
    )
    """
    Compute linear speed of points on rotating cylinder
    θ = atan(y/x)
    ẋ = -R*θ̇*sin(θ)
    ẏ =  R*θ̇*cos(θ)

    Note: this is a specialized case for development only
    """
    xb = body.xb
    ub = body.ub
    motion = body.motion
    R = sqrt(xb[1, 1]^2 + xb[1, 2]^2)
    println(R)
    θ = atan.(xb[:, 2], xb[:, 1])
    ub[:, 1] .= -R*motion.θ̇*sin.(θ)
    ub[:, 2] .=  R*motion.θ̇*cos.(θ)
    return ub
end


function move_body!(
    ::Type{MotionFunction},
    body::RigidBody,
    t::Float64
    )
    se2 = get_transformation(body.motion, t)
    #println(ub[1, 1])
    #println(se2.Ru[1, :] * ub[1, :])
    for i=1:size(body.xb, 1)
        @views body.xb[i, :] .= se2.pos .+ se2.Rx*body.x0[i, :]
        @views body.ub[i, :] .= se2.vel .+ se2.Ru*body.x0[i, :]
    end
    #println(se2.Ru)
    #return ub
    #xb .= se2.pos .+ se2.Rx*xb
    #ub .= se2.vel .+ se2.Ru*ub
end
