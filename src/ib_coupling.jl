function delta_h( rf, rb , dr )
    """
    %take points on the flow domain (r) that are within the support
    %(supp) of the IB points (rb), and evaluate
    %               delta( abs(r - rb) )

    %Currently uses the Yang3 smooth delta function (see Yang et al, JCP,
    %2009), which has a support of 6*h (3*h on each side)
    """

    # Note: the result is delta * h
    r = abs.( rf .- rb );
    r1 = r./dr;
    r2 = r1.*r1;
    r3 = r2.*r1;
    r4 = r3.*r1;

    del_h = zeros(size(r1));
    a5 = zeros(size(r1))
    a8 = zeros(size(r1))
    a6 = zeros(size(r1))
    a9 = zeros(size(r1))
    a1 = zeros(size(r1))
    a7 = zeros(size(r1))

    ind = r1 .< 1.0  # Have to avoid domain errors here
    a5[ind] .= @. asin((1.0/2.0)*sqrt(3.0)*(2.0*r1[ind] - 1.0));
    a8[ind] .= @. sqrt(1.0 - 12.0*r2[ind] + 12.0*r1[ind]);

    del_h .+= @. ind.*(
        4.166666667e-2*r4 + (-0.1388888889+3.472222222e-2*a8)*r3 +
        (-7.121664902e-2 - 5.208333333e-2*a8 + 0.2405626122*a5)*r2 +
        (-.2405626122*a5 - 0.3792313933 + 0.1012731481*a8)*r1 + 8.0187537413e-2*a5 -
        4.195601852e-2*a8 + 0.6485698427 )

    ind = (r1.>1.0).*(r1.<=2.0)
    a6[ind] .= @. asin((1.0/2.0)*sqrt(3.0)*(-3.0+2.0*r1[ind]));
    a9[ind] .= @. sqrt(-23.0+36.0*r1[ind]-12.0*r2[ind]);

    del_h .+= @. ind.*(
        -6.250000000e-2*r4 + (0.4861111111 - 1.736111111e-2*a9).*r3 +
        (-1.143175026+7.812500000e-2*a9 - 0.1202813061*a6)*r2 +
        (0.8751991178+0.3608439183*a6 -0.1548032407*a9)*r1 - 0.2806563809*a6 +
        8.22848104e-3 + 0.1150173611*a9 )

    ind = (r1.>2.0).*(r1.<=3.0);
    a1[ind] .= @. asin((1.0/2.0*(2.0*r1[ind]-5.0))*sqrt(3.0))
    a7[ind] .= @. sqrt(-71.0-12.0*r2[ind]+60.0*r1[ind])

    del_h .+= @. ind.*(
        2.083333333e-2*r4 + (3.472222222e-3*a7 -0.2638888889)*r3 +
        (1.214391675 - 2.604166667e-2*a7 + 2.405626122e-2*a1)*r2 +
        (-0.1202813061*a1 - 2.449273192 + 7.262731481e-2*a7)*r1 + 0.1523563211*a1 +
        1.843201677 - 7.306134259e-2*a7 )

    return del_h
end



function coupling_mat( grid::T, bodies::Array{<:Body, 1}; supp=6.0 ) where T <: Grid
    """
    Return the regularization matrix E' (scaled to be the transpose of E),
    which takes quantities on the IB and smears them to the flow domain.

    NOTE: E should only be precomputed for Static motions!

    E is the interpolation matrix, which goes from the flow domain to the IB
        This depends on the position of the body points, but not their velocities

    200 x 200 grid
    Original:
        99.829 ms (19585 allocations: 116.02 MiB)
    Streamlined:
        26.615 ms (15461 allocations: 97.54 MiB)
    Direct CSC construction:
        15.914 ms (15902 allocations: 29.88 MiB)
    """
    m = grid.nx;
    n = grid.ny;
    len = grid.len;
    offx = grid.offx;
    offy = grid.offy;

    del = grid.h;  # Size of uniform grid cell

    # Get size of ET
    # Also need later: y index starts after all the x-vels
    n_add = get_velx_ind( m-1, n, grid )
    nrows = n_add + get_vely_ind( m, n-1, grid );
    nb, nf = get_body_info(bodies)
    ncols = sum(nf);

    tally = 0; # Used to keep a tally of the column we're on

    # for x-vels: x and y points on physical grid (for fluid domain)
    xu = (del : del : (m-1) * del) .- offx;
    yu = (del/2 : del : (n-1/2) * del ) .- offy;

    # for y-vels: points on physical grid (for fluid domain)
    xv = (del/2. : del : (m-1/2) * del ) .- offx;
    yv = (del : del : (n-1)*del ) .- offy;

    colptr = ones(Int64, ncols+1);
    rowval = Array{Int64}(undef, 0, 1);
    nzval = Array{Float64}(undef, 0, 1);

    for i = 1:length(bodies)

        # x and y points on physical grid (for IB)
        xb_x = bodies[i].xb[:, 1];
        xb_y = bodies[i].xb[:, 2];

        # For each IB point, add nonzero weights...
        #--rows corresponding to x-vels
        for j = 1:nb[i]
            # Find x-points of fluid domain within support
            ind_x =  findall( abs.( xu .- xb_x[j] ) .<= supp * del )

            # Find y-points within support
            ind_y = findall( abs.( yu .- xb_y[j] ) .<= supp * del )

            # rows of ET to add to (first grid level)
            velx_ind = (m-1).*(ind_y' .- 1) .+ ind_x

            # entries to put into columns
            del_h = delta_h( xu[ind_x], xb_x[j], del) .*
                    delta_h( yu[ind_y], xb_y[j], del)'

            colptr[tally+j+1] = colptr[tally+j] + length(velx_ind)
            rowval = [rowval; velx_ind[:]]
            nzval = [nzval; del_h[:]]

            # E[tally+j, velx_ind[:]] .+= del_h[:];  # without direct CSC construction
        end

        tally += nb[i];

        #--rows corresponding to y-vels
        for j = 1:nb[i]

            # Find x-points of fluid domain within support
            ind_x =  findall( abs.( xv .- xb_x[j] ) .<= supp * del )

            # Find y-points within support
            ind_y = findall( abs.( yv .- xb_y[j] ) .<= supp * del )

            # rows of ET to add to
            vely_ind = n_add .+ m .* (ind_y' .- 1) .+ ind_x

            # entries to put into columns
            del_h = delta_h( xv[ind_x], xb_x[j], del) .*
                    delta_h( yv[ind_y], xb_y[j], del)'

            colptr[tally+j+1] = colptr[tally+j] + length(vely_ind)
            rowval = [rowval; vely_ind[:]]
            nzval = [nzval; del_h[:]]

            # E[tally+j, vely_ind[:]] .+= del_h[:];  # without direct CSC construction
        end

        #Keep a tally of columns we've used to now
        tally += nb[i];
    end

    return SparseMatrixCSC(nrows, ncols, colptr, rowval[:], nzval[:])'
end
