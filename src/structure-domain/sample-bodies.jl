"""
Support for constructing basic bodies:
    * Flat 1D plate
    * Circular cylinder
    * 4-digit NACA airfoils
"""

"""
    make_plate(L, α, h, x0, y0; motion=Static(), n=0)

Build plate of length L at AoA α based on grid spacing of h.

If n = 0, choose so that ds=2h.  Also note that α here is different from
the free-stream flux definition, i.e. you can either incline the plate or change
the angle of the free-stream flow (but you probably don't want to do both)
"""
function make_plate(L, α, h, x0, y0; motion=Static(), n=0)
    # Get # of points such that ds = 2h
    if (n==0)
        n = Int( floor( L / h / 2 ) ) + 1;
    end

    spt = L.*(0:(n-1))/(n-1);  # Range (0, L)
    xhat = spt*cos.(-α);
    yhat = spt*sin.(-α);

    xb = [xhat.+x0  yhat.+y0];

    # sanity check: make sure ds is equal to 2 * h
    ds = sqrt( (xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2 ) ;

    return RigidBody(motion, xb, copy(xb), 0.0*xb, fill(ds, n))
end

"""
    make_cylinder(r, h, x0, y0; motion=Static(), n=0)

Build cylinder of radius r based on grid spacing of h.

If n \neq 0, choose so that ds=2h.
"""
function make_cylinder(r, h, x0, y0; motion=Static(), n=0)
    circum = 2 * π * r; #  Circumference of the circle

    if ismissing(h)==true
        h = 2.0/Re #default to a grid Re of 2
    end

    # Get # of points such that ds = 2h
    if (n==0)
        n = Int( floor( circum / h / 2 ) );
    end

    int =  2*π/n ;
    spt = 0:int:(n-1)*int;
    xhat = r.*cos.(spt);
    yhat = r.*sin.(spt);

    xb = [xhat.-x0  yhat.-y0];

    # sanity check: make sure ds is equal to 2 * h
    ds = sqrt( (xhat[2] - xhat[1])^2 + (yhat[2] - yhat[1])^2 ) ;

    return RigidBody(motion, xb, copy(xb), 0.0*xb, fill(ds, n))
end

"""
    make_naca(x0, N, spec; motion=Static())

Generate 4-digit NACA airfoil based on string `spec`.
"""
function make_naca(x0, N, spec; motion=Static())

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

"""
    naca(x, spec)

Compute points on 4-digit NACA airfoils

x - x/c, so that x ∈ (0, 1)
"""
function naca(x, spec)
    # First, break down spec
    m = parse(Int, spec[1])/100.    # Maximum camber
    p = parse(Int, spec[2])/10.     # Location of max. camber
    t = parse(Int, spec[3:4])/100.  # Maximum thickness

    yc = zeros(size(x));
    if (m > 0) && (p > 0)  # Cambered airfoil
        max_idx = findall(x.<p)[end]

        # Compute mean camber line
        yc[1:max_idx] .= (m/p^2).*(2*p.*x[1:max_idx] .- x[1:max_idx].^2)
        yc[max_idx+1:end] .= (m/(1-p)^2).*( 1-2*p .+ 2*p.*x[max_idx+1:end] .- x[max_idx+1:end].^2)

        dyc_dx = zeros(size(x));
        dyc_dx[1:max_idx] .= (2*m/p^2).*(p .- x[1:max_idx]);
        dyc_dx[max_idx+1:end] .= (2*m/(1-p^2)).*(p .- x[max_idx+1:end]);
        θ = atan.(dyc_dx)
    else
        θ = zeros(size(x))
    end

    yt = 5*t*(0.2969*sqrt.(x) - 0.1260.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1036.*x.^4 );

    xU = x .- yt.*sin.(θ)
    xL = x .+ yt.*sin.(θ)
    yU = yc .+ yt.*cos.(θ)
    yL = yc .- yt.*cos.(θ)

    return yc, xU, xL, yU, yL
end

"""
    sym_naca(x, spec)

Function naca(x, t) for points on symmetric 4-digit airfoils
"""
function sym_naca(x, spec)
    sym_naca(x, t) = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1036*x^4 );
end
