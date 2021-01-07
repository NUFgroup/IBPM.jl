


abstract type ExplicitScheme end


struct AdamsBashforth <: ExplicitScheme
    dt::Float64
    β::Array{Float64, 1}
end

function AB2(dt::Float64)
    """
    Initialize second-order Adams-Bashforth scheme
    """
    return AdamsBashforth(dt, [1.5, -0.5])
end
