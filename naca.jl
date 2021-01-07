
function naca(x, spec)
    """
    Compute points on 4-digit NACA airfoils
    x - x/c, so that x ∈ (0, 1)
    """
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

    #return xU, xL, yU, yL
    return yc, xU, xL, yU, yL
end

function sym_naca(x, spec)
    """
    Compute points on symmetric 4-digit airfoils
    """
    sym_naca(x, t) = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1036*x^4 );
end
