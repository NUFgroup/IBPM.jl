"""discrete delta function used to relate flow to structure quantities"""

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

"""
Use the discrete delta function to build ET, which regularizes info from the
    immersed surface onto the flow domain
"""
# TODO: optimize or (better yet) make an implicit operator
function reg_mats( grid::T, bodies::Array{<:Body, 1}; supp=6.0 ) where T <: Grid
    """
    Return the regularization matrix ET (scaled to be the transpose of E),
    which takes quantities on the IB and smears them to the flow domain.

    This depends on the position of the body points, but not their velocities
    """
    m = grid.nx;
    n = grid.ny;
    len = grid.len;
    offx = grid.offx;
    offy = grid.offy;

    del = grid.h;  # Size of uniform grid cell

    # Get size of ET
    nrows = get_velx_ind( m-1, n, grid ) +
        get_vely_ind( m, n-1, grid );
    nb, nf = get_body_info(bodies)
    ncols = sum(nf);

    ET = spzeros( nrows, ncols );
    ncol_tally = 0; # Used to keep a tally of the column we're on

    for i = 1:length(bodies)

        xb = bodies[i].xb;
        #--rows corresponding to x-vels

            #x and y points on physical grid (for fluid domain)
            xu = (del : del : (m-1)*del ) .- offx;
            yu = (del/2 : del : (n-1/2) * del ) .- offy;

            # x and y points on physical grid (for IB)
            xb_x = xb[:, 1];
            xb_y = xb[:, 2];

            # find x points within support
            xu_r = repeat( xu, length(xb_x) );
            xb_rx = repelem( xb_x, length(xu) )
            ind_supp_x = ( abs.( xu_r .- xb_rx ) .<= supp * del )

            # Get index of these points for fluid grid
            i_supp = Int.( round.( (xu_r[ind_supp_x] .+ offx) / del ) )

            # Get index of these points for IB:
            i_supp_xbx = Int.( round.( ind_supp_x .* repelem( 1:nb[i], length(xu) ) ) );
            i_supp_xbx = i_supp_xbx[ i_supp_xbx .!= 0 ];

            # find y points within support
            yu_r = repeat( yu, length(xb_y));
            xb_ry = repelem( xb_y, length(yu) );
            ind_supp_y = ( abs.( yu_r .- xb_ry ) .<= supp * del )

            # Get y-index of these points for fluid grid
            j_supp = Int.( round.( (yu_r[ind_supp_y] .+ offy) / del .+ 0.5  ) )

            # Get index of these points for IB:
            j_supp_xby = Int.( round.( ind_supp_y .* repelem( 1:nb[i], length(yu) ) ) )
            j_supp_xby = j_supp_xby[ j_supp_xby .!= 0 ]

            # For each IB point, add nonzero weights...
            for j = 1:nb[i]

                # x-indices on IB corresponding to current body point
                ind_xbx = findall(i_supp_xbx .== j)

                # x-indices on flow grid that are within support of IB point
                ind_x = i_supp[ ind_xbx ];

                # y-indices on IB corresponding to current body point
                ind_xby = findall(j_supp_xby .== j);

                # y-indices on flow grid that are within support of IB point
                ind_y = j_supp[ ind_xby ];

                # Combine flow indices
                indvelx = repeat( ind_x, length(ind_y));
                indvely = repelem( ind_y, length(ind_x) );

                # rows of ET to add to
                velx_ind = get_velx_ind( indvelx, indvely, grid );

                # columns to add to
                xb_ind = j .* ones(Int32, size(velx_ind ) );

                # entries to put into columns
                del_h = delta_h( xu[indvelx], xb_x[xb_ind], del) .*
                    delta_h( yu[indvely], xb_y[xb_ind], del)


                # Add to ET:
                # y index starts after all the x-vels
                ET .+= sparse( velx_ind, ncol_tally .+ xb_ind,
                               del_h, nrows, ncols );
            end

        #--rows corresponding to y-vels

            # x and y points on physical grid (for fluid domain)
            xv = (del/2. : del : (m-1/2) * del ) .- offx;
            yv = (del : del : (n-1)*del ) .- offy;

            # x and y points on physical grid (for IB)
            xb_x = xb[ 1 : nb[i] ];
            xb_y = xb[ 1 + nb[i] : 2*nb[i] ];

            # find x points within support
            xv_r = repeat( xv, length(xb_x));
            xb_rx = repelem( xb_x, length(xv) );
            ind_supp_x = ( abs.( xv_r .- xb_rx ) .<= supp * del )

            # Get index of these points for fluid grid
            i_supp = Int.( round.( (xv_r[ind_supp_x] .+ offx)/del .+ 0.5 ) )

            # Get index of these points for IB:
            i_supp_xbx = Int.( round.( ind_supp_x .* repelem( 1:nb[i], length(xv) ) ) )
            i_supp_xbx = i_supp_xbx[ i_supp_xbx .!= 0 ]

            # find y points within support
            yv_r = repeat( yv, length(xb_y))
            xb_ry = repelem( xb_y, length(yv) )
            ind_supp_y = ( abs.( yv_r .- xb_ry ) .<= supp * del )

            # Get y-index of these points for fluid grid
            j_supp = Int.( round.( (yv_r[ind_supp_y] .+ offy) / del ) )

            # Get index of these points for IB:
            j_supp_xby = Int.( round.( ind_supp_y .* repelem( 1:nb[i], length(yv) ) ) )
            j_supp_xby = j_supp_xby[ j_supp_xby .!= 0 ]

            #println(sum(j_supp))
            # For each IB point, add nonzero weights...
            for j = 1:nb[i]

                # x-indices on IB corresponding to current body point
                ind_xbx = findall(i_supp_xbx .== j);

                # x-indices on flow grid that are within support of IB point
                ind_x = i_supp[ ind_xbx ];

                # y-indices on IB corresponding to current body point
                ind_xby = findall(j_supp_xby .== j);

                # y-indices on flow grid that are within support of IB point
                ind_y = j_supp[ ind_xby ];

                # Combine flow indices
                indvelx = repeat( ind_x, length(ind_y) );
                indvely = repelem( ind_y, length(ind_x) );

                # rows of ET to add to
                vely_ind = get_vely_ind( indvelx, indvely, grid );

                # columns to add to
                xb_ind = j .* ones(Int32, size(vely_ind ) );

                # entries to put into columns
                del_h = delta_h( xv[indvelx], xb_x[xb_ind], del) .*
                    delta_h( yv[indvely], xb_y[xb_ind], del)

                # Add to ET:
                # y index starts after all the x-vels
                n_add = get_velx_ind( m-1, n, grid )
                ET .+= sparse( n_add .+ vely_ind, xb_ind .+ ncol_tally .+ nb[i],
                    del_h, nrows, ncols );
            end


        #Keep a tally of columns we've used to now
        ncol_tally = ncol_tally + length(xb);

    end
    return ET
end

"""
Compute the action of the matrix, B, that is used to solve for the surface
stresses that lead to a velocity that satisfies the no-slip BCs
"""
function b_times!(x::Array{Float64, 2},
                  z::Array{Float64, 2},
                  Ainv::LinearMap,
                  model::IBModel{UniformGrid, RigidBody{Static}},
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.

    % (B arises from an LU factorization of the full system)

    Note ψ is just a dummy work array for circ2_st_vflx
    """
    # -- get circulation from surface stress  circ = Ainv * R * E' * z
    #     We don't include BCs for Ainv because ET*z is compact
    #circ = Ainv( mats.R*(mats.ET*z), dt, model );
    mul!(Γ, Ainv, model.mats.RET*z)

    #-- get vel flux from circulation
    #vflx, _ = circ2_st_vflx( circ, model );
    circ2_st_vflx!( ψ, q, Γ, model );

    #--Interpolate onto the body and scale by h
    #x = (model.mats.E*vflx) / model.grid.h;
    mul!(x, model.mats.E, q)
    rmul!(x, 1/model.grid.h)
end



# 44.733 ms (321 allocations: 48.80 MiB) for one evaluation
function b_times!(x::Array{Float64, 2},
                  z::Array{Float64, 2},
                  Ainv::LinearMap,
                  model::IBModel{MultiGrid, RigidBody{Static}},
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.

    % (B arises from an LU factorization of the full system)

    Note ψ is just a dummy variable for computing velocity flux
        Also this only uses Ainv on the first level
    """

    # --Initialize
    grid = model.grid
    mats = model.mats
    m = grid.nx;
    n = grid.ny;
    nΓ = grid.nΓ   # Number of circulation points
    mg = grid.mg

    # -- get circulation from surface stress

    # Get circ on 1st grid level
    # TODO: in place with @view macro
    Γ[:, 1] = Array( Ainv * (mats.RET*z) );
    # We don't include BCs from coarse grid for Ainv because ET*z is compact

    # Coarsify circulation to second grid level to get BCs for stfn
    #Γ[:, 2] = coarsify( Γ[:,1], Γ[:,2], grid );
    @views coarsify!( Γ[:,1], Γ[:,2], grid );

    #-- get vel flux from circulation

    # only need to work with 2 grid levels here
    #vflx, _ = circ2_st_vflx( circ, 2, model);
    circ2_st_vflx!( ψ, q, Γ, model, 2 );

    #--Interpolate onto the body and scale by h
    mul!(x, mats.E, q[:, 1])
    rmul!(x, 1/grid.h)

end
