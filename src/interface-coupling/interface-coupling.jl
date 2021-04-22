"""
    δh( rf, rb , dr )

Discrete delta function used to relate flow to structure quantities
"""
function δh( rf, rb , dr )
    """
    %take points on the flow domain (r) that are within the support
    %(supp) of the IB points (rb), and evaluate
    %               delta( abs(r - rb) )

    %Currently uses the Yang3 smooth delta function (see Yang et al, JCP,
    %2009), which has a support of 6*h (3*h on each side)

    Note that this gives slightly different answers than Fortran at around 1e-4,
        apparently due to slight differences in the floating point arithmetic.
        As far as I can tell, this is what sets the bound on agreement between
        the two implementations.  It's possible this might be improved with
        arbitrary precision arithmetic (i.e. BigFloats), but at least it
        doesn't seem to be a bug.
    """

    # Note: the result is delta * h
    r = abs( rf - rb );
    r1 = r/dr; r2 = r1*r1; r3 = r2*r1; r4 = r3*r1;

    if (r1 <= 1.0)
        #println("r < 1")
        a5 = asin((1.0/2.0)*sqrt(3.0)*(2.0*r1 - 1.0));
        a8 = sqrt(1.0 - 12.0*r2 + 12.0*r1);

        del_h = 4.166666667e-2*r4 + (-0.1388888889+3.472222222e-2*a8)*r3 +
            (-7.121664902e-2 - 5.208333333e-2*a8 + 0.2405626122*a5)*r2 +
            (-.2405626122*a5 - 0.3792313933 + 0.1012731481*a8)*r1 + 8.0187537413e-2*a5 -
            4.195601852e-2*a8 + 0.6485698427

    elseif (r1 <= 2.0)
        #println("r < 2")
        a6 = asin((1.0/2.0)*sqrt(3.0)*(-3.0+2.0*r1));
        a9 = sqrt(-23.0+36.0*r1-12.0*r2);

        del_h = -6.250000000e-2*r4 + (0.4861111111 - 1.736111111e-2*a9).*r3 +
            (-1.143175026+7.812500000e-2*a9 - 0.1202813061*a6)*r2 +
            (0.8751991178+0.3608439183*a6 -0.1548032407*a9)*r1 - 0.2806563809*a6 +
            8.22848104e-3 + 0.1150173611*a9

    elseif (r1 <= 3.0 );
        #println("r < 3")
        a1 = asin((1.0/2.0*(2.0*r1-5.0))*sqrt(3.0))
        a7 = sqrt(-71.0-12.0*r2+60.0*r1)


        del_h = 2.083333333e-2*r4 + (3.472222222e-3*a7 -0.2638888889)*r3 +
        (1.214391675 - 2.604166667e-2*a7 + 2.405626122e-2*a1)*r2 +
        (-0.1202813061*a1 - 2.449273192 + 7.262731481e-2*a7)*r1 + 0.1523563211*a1 +
        1.843201677 - 7.306134259e-2*a7
    else
        del_h = 0.0
    end

    return del_h
end

"""
Return the interpolation matrix E, which interpolates flow quantities
to the immersed boundary.  This depends on the position of the body points,
but not their velocities

The transpose of E is the regularization matrix which smears quantities
from the IB to the flow domain.
"""
function setup_reg( grid::T, bodies::Array{<:Body, 1}; supp=6 ) where T <: Grid
    nx = grid.nx; ny = grid.ny; h = grid.h;  # Size of uniform grid cell

    # Stack all body points together... don't need to distinguish them here
    xb = vcat( [body.xb for body in bodies]...)
    nb = size(xb, 1)
    nf = 2*nb

    supp_idx = -supp:supp;
    #l=supp_idx.+(supp+1); m=supp_idx'.+(supp+1)  # Shift indices to index into "weight" at 1
    weight = zeros( nf, 2, 2*supp+1, 2*supp+1 )

    # Nearest indices of body relative to grid
    body_idx = zeros(Int, size(xb))
    body_idx[:, 1] .= @. Int(floor( (xb[:, 1]+grid.offx)/h ))
    body_idx[:, 2] .= @. Int(floor( (xb[:, 2]+grid.offy)/h ))

    k, l, m = 1, 6, 5

    x = @. h*(body_idx[k, 1]-1+supp_idx)-grid.offx       # grid location x
    y = @. h*(body_idx[k, 2]-1+supp_idx')-grid.offy       # grid location y

    # get regularized weight near IB points (u-vel points)
    for k=1:nb
        x = @. h*(body_idx[k, 1]-1+supp_idx)-grid.offx       # grid location x
        y = @. h*(body_idx[k, 2]-1+supp_idx')-grid.offy       # grid location y

        @. weight[k, 1, :, :] = δh(x, xb[k, 1], h) * δh(y+h/2, xb[k, 2], h)
        @. weight[k, 2, :, :] = δh(x+h/2, xb[k, 1], h) * δh(y, xb[k, 2], h)
    end

    " Matrix E' "
    function reg!(q, fb)
        q .*= 0.0
        fb = reshape(fb, nb, 2)
        qx, qy = grid.split_flux(q)
        for k=1:nb
            i=body_idx[k, 1].+supp_idx; j=body_idx[k, 2].+supp_idx
            @views qx[i, j] += weight[k, 1, :, :]*fb[k, 1]
            @views qy[i, j] += weight[k, 2, :, :]*fb[k, 2]
        end
        fb = reshape(fb, 2*nb, 1)
        return nothing
    end

    " Matrix E "
    function regT!(fb, q)
        fb .*= 0.0
        fb = reshape(fb, nb, 2)
        qx, qy = grid.split_flux(q)
        for k=1:nb
            i=body_idx[k, 1].+supp_idx; j=body_idx[k, 2].+supp_idx
            fb[k, 1] += sum( weight[k, 1, :, :].*qx[ i, j ] )
            fb[k, 2] += sum( weight[k, 2, :, :].*qy[ i, j ] )
        end
        fb = reshape(fb, 2*nb, 1)
        return nothing
    end

    E = LinearMap( regT!, reg!, nf, grid.nq; ismutating=true  )
    return E
end
