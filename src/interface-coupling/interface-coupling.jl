"""discrete delta function used to relate flow to structure quantities"""
function delta_h( rf, rb , dr; tol=1e-12 )
    """
    %take points on the flow domain (r) that are within the support
    %(supp) of the IB points (rb), and evaluate
    %               delta( abs(r - rb) )

    %Currently uses the Yang3 smooth delta function (see Yang et al, JCP,
    %2009), which has a support of 6*h (3*h on each side)
    """

    # Note: the result is delta * h
    r = abs( rf - rb );
    r1 = r/dr;
    r2 = r1*r1;
    r3 = r2*r1;
    r4 = r3*r1;

    if (r1 <= 1.0)
        a5 = asin((1.0/2.0)*sqrt(3.0)*(2.0*r1 - 1.0));
        a8 = sqrt(1.0 - 12.0*r2 + 12.0*r1);

        del_h = 4.166666667e-2*r4 + (-0.1388888889+3.472222222e-2*a8)*r3 +
            (-7.121664902e-2 - 5.208333333e-2*a8 + 0.2405626122*a5)*r2 +
            (-.2405626122*a5 - 0.3792313933 + 0.1012731481*a8)*r1 + 8.0187537413e-2*a5 -
            4.195601852e-2*a8 + 0.6485698427

    elseif (r1 <= 2.0)
        a6 = asin((1.0/2.0)*sqrt(3.0)*(-3.0+2.0*r1));
        a9 = sqrt(-23.0+36.0*r1-12.0*r2);

        del_h = -6.250000000e-2*r4 + (0.4861111111 - 1.736111111e-2*a9).*r3 +
            (-1.143175026+7.812500000e-2*a9 - 0.1202813061*a6)*r2 +
            (0.8751991178+0.3608439183*a6 -0.1548032407*a9)*r1 - 0.2806563809*a6 +
            8.22848104e-3 + 0.1150173611*a9

    elseif (r1 <= 3.0 );
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
    nx = grid.nx; ny = grid.ny;
    len = grid.len;
    offx = grid.offx;
    offy = grid.offy;

    del = grid.h;  # Size of uniform grid cell
    d2 = del/2

    # Stack all body points together... don't need to distinguish them here
    xb = vcat( [body.xb for body in bodies]...)
    nb = size(xb, 1)
    nf = 2*nb

    weight = zeros( nf, (2*supp + 1)^2 )
    indexx = zeros(Int64, nb, 2 )

    # Nearest indices of body relative to grid
    for i=1:nb
        indexx[i, 1] = Int(round( (xb[i, 1]+grid.offx)/del ))
        indexx[i, 2] = Int(round( (xb[i, 2]+grid.offy)/del ))
    end

    # get regularized weight near ib points (u-vel points)
    for i=1:nb
        next=0
        for l=-supp:supp
            for k=-supp:supp
                x = del*(indexx[i, 1]-1+k)-grid.offx       # grid location x
                y = del*(indexx[i, 2]-1+l)-grid.offy+d2 # grid location y
                next += 1
                weight[i, next] = delta_h(x, xb[i, 1], del) *
                                 delta_h(y, xb[i, 2], del)

                x = del*(indexx[i, 1]-1+k)-grid.offx+d2  # grid location x
                y = del*(indexx[i, 2]-1+l)-grid.offy # grid location y
                weight[i+nb,next] = delta_h(x, xb[i, 1], del) *
                                    delta_h(y, xb[i, 2], del)
            end
        end
    end

    function reg!(q, fb)
        #q = zeros(grid.nq)
        q .*= 0.0
        for k=1:nb
            i = indexx[k, 1]
            j = indexx[k, 2]
            next = 0
            for l=-supp:supp
                for p=-supp:supp
                    next += 1
                    q[grid.u_idx(i+p,j+l)] += weight[k,next]*fb[k]
                    q[grid.v_idx(i+p,j+l)] += weight[k+nb,next]*fb[k+nb]
                end
            end
        end
    end

    function regT!(fb, q)
        fb .*= 0.0
        for k=1:nb
            i = indexx[k]
            j = indexx[k+nb]
            next = 0
            for l=-supp:supp
                for p=-supp:supp
                    next += 1
                    fb[k] += weight[k,next]*q[ grid.u_idx(i+p,j+l) ]
                    fb[k+nb] += weight[k+nb,next]*q[ grid.v_idx(i+p,j+l) ]
                end
            end
        end
    end

    E = LinearMap( regT!, reg!, nf, grid.nq; ismutating=true  )
    return E
end
