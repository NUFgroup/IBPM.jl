using FFTW
using LinearAlgebra

# Shorthand for discrete sine transform
function dst(b)
    #return FFTW.r2r(b, FFTW.RODFT10, 2)
    return 0.25*FFTW.r2r(b, FFTW.RODFT00, [1, 2])
end

function get_dst_plan(b)
    # Generate plans needed for optimized DST inversion
    p = 0.25*FFTW.plan_r2r(b, FFTW.RODFT00, [1, 2], flags=FFTW.UNALIGNED)

    # Preallocated work array
    w = zeros(size(b))

    return p, w
end

# Original:    8.409 ms (264 allocations: 11.85 MiB)
# 2D:          5.323 ms (128 allocations: 4.56 MiB)
function dst_inv(b, Λ)
    return dst( dst(b) ./Λ' )
    #return dst( transpose( dst( dst( transpose( dst( b ) ) ) ./ Λ ) ) )
end

# Original:    5.083 ms (16 allocations: 6.37 MiB)
# 2D:          4.312 ms (6 allocations: 2.73 MiB)
function dst_inv(b, Λ, dst_plan)
    # p - plan for dst of b
    # pT - plan for dst of p(b)'
    #p, pT = dst_plan[1:2];
    #return p * transpose( pT * ( ( pT * transpose(p * b) ) ./ Λ ) )

    p = dst_plan[1];
    return p * ( (p * b) ./ Λ' )
end

function dst_inv!(x, b, Λ, dst_plan)
    """
    Optimized DST inversion with pre-allocated arrays
    """
    # p - plan for dst of b
    # w - work array (preallocated memory)
    p, w = dst_plan
    mul!(w, p, b)  # dst(b)
    w ./= transpose(Λ);
    mul!(x, p, w)  # dst( (dst(b) ./ Λ') )
end
