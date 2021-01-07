"""
Functions to create sparse matrices for IBPM

Note: these matrices are "static" in the sense that they act on the field
    but are not used in the time-stepping scheme (i.e. they don't depend on Δt)
    Those matrices are in 'timestepper_utils.jl'
"""

# Indexing helper functions
x1_ind(m::Int, n::Int) = repeat(1:m-1, n-1)
x2_ind(m::Int, n::Int) = repeat(2:m, n-1)
y1_ind(m::Int, n::Int) = repelem(1:n-1, m-1)   # repelem is a clone of MATLAB's function.. in ib_matutils.jl
y2_ind(m::Int, n::Int) = repelem(2:n, m-1)

function get_C( grid::T ) where T <: Grid
    m = grid.nx;
    n = grid.ny;

    nrows = get_velx_ind( m-1, n, grid ) +
            get_vely_ind( m, n-1, grid ) ;
    ncols = get_vort_ind( m-1, n-1, grid ) ;
    C = spzeros(nrows, ncols)

    # --First build block corresponding to x-velocities

    # vorticity points above x-velocity point
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);

    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    C .+= sparse( velx_ind, vort_ind,
              ones(size(vort_ind)), nrows, ncols);

    #vorticity points below x-velocity point (only the velx_ind chages)
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);

    C .-= sparse( velx_ind, vort_ind,
              ones(size(vort_ind)), nrows, ncols);

    #--Now build y-velocity block

    # rows start at end of x-velocity block:
    n_add = get_velx_ind( m-1, n, grid );

    # vorticity points to the right of y-velocity point
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);

    C .-= sparse(n_add .+ vely_ind, vort_ind,
                ones(size(vort_ind)), nrows, ncols);

    # vorticity points to the left of y-velocity point
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);

    C .+= sparse(n_add .+ vely_ind, vort_ind,
                ones(size(vort_ind)), nrows, ncols);

    return C
end



function get_Q( grid::T ) where T <: Grid
    """
    %Build averaging operator Q = [Qx; Qy]. The x-velocity block (Qx) takes y
    %velocities and averages them onto the x-velocity edges, and the y-velocity
    %block (Qy) takes x velocities and averages them onto y-velocity edges.

    WHAT DOES THIS MEAN???
    %Note: this does not give the Q used in the code. The Q used in the code is
    %postmultiplies Minv to the Q obtained here.

    %Inputs: grid_parms -- data structure containing m (number of points in x
    %dirn), n (number of points in y dirn), and mg (number of grid levels)

    """

    m = grid.nx;
    n = grid.ny;

    nrows = get_velx_ind( m-1, n, grid ) +
            get_vely_ind( m, n-1, grid ) ;
    ncols = nrows ;
    Q = spzeros(nrows, ncols)  # Initialize


    # y-velocity index starts at the end of all x-velocities
    n_add = get_velx_ind( m-1, n, grid );

    #--Block corresponding to x-velocities

    # y-vels to bottom left of current x-vel point
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);

    Q .-= 0.25*sparse( velx_ind, n_add .+ vely_ind,
            ones(size(vely_ind)), nrows, ncols);

    # y-vels to bottom right of current x-vel point
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);

    Q .-= 0.25*sparse( velx_ind, n_add .+ vely_ind,
            ones(size(vely_ind)), nrows, ncols);

    # y-vels to top left of current x-vel point
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);

    Q .-= 0.25*sparse( velx_ind, n_add .+ vely_ind,
            ones(size(vely_ind)), nrows, ncols);

    #y-vels to top right of current x-vel point
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);

    Q .-= 0.25*sparse( velx_ind, n_add .+ vely_ind,
            ones(size(vely_ind)), nrows, ncols);


    #--Block corresponding to y-velocities

    # x-vels to bottom left of current y-vel point
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);

    Q .+= 0.25*sparse( n_add .+ vely_ind, velx_ind,
        ones(size(vely_ind)), nrows, ncols);

    # x-vels to bottom right of current y-vel point
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);

    Q .+= 0.25*sparse( n_add .+ vely_ind, velx_ind,
        ones(size(vely_ind)), nrows, ncols);

    # x-vels to top left of current y-vel point\
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);

    Q .+= 0.25*sparse( n_add .+ vely_ind, velx_ind,
        ones(size(vely_ind)), nrows, ncols);

    # x-vels to top right of current y-vel point
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);

    Q .+= 0.25*sparse( n_add .+ vely_ind, velx_ind,
        ones(size(vely_ind)), nrows, ncols);

    return Q
end




function get_W( grid::T ) where T <: Grid
    """
    %Build averaging operator W = [Wx; Wy]. The x-velocity block (Wx) takes
    %vorticity and averages it onto the x-velocity edges, and the y-velocity
    %block (Wy) averages it onto y-velocity edges.

    %Note that W is not the full W used in the code. The W in the code is M *
    %the W computed here.

    %Inputs: grid_parms -- data structure containing m (number of points in x
    %dirn), n (number of points in y dirn), and mg (number of grid levels)
    """

    m = grid.nx;
    n = grid.ny;

    nrows = get_velx_ind( m-1, n, grid ) +
            get_vely_ind( m, n-1, grid ) ;
    ncols = get_vort_ind( m-1, n-1, grid ) ;
    W = spzeros(nrows, ncols)  # Initialize

    # y-velocity index starts at the end of all x-velocities
    n_add = get_velx_ind( m-1, n, grid );

    #-- x-vel block
    # vorticity points above x-velocity point
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);
    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    W .+= 0.5*sparse( velx_ind, vort_ind,
        ones(size(vort_ind)), nrows, ncols);

    # vorticity points below x-velocity point
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);
    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    W .+= 0.5*sparse( velx_ind, vort_ind,
        ones(size(vort_ind)), nrows, ncols);


    #-- y-vel block
    #vorticity points to the right of y-velocity point
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);
    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    W .+= 0.5*sparse( n_add .+ vely_ind, vort_ind,
        ones(size(vort_ind)), nrows, ncols);

    #vorticity points to the left of y-velocity point
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);
    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    W .+= 0.5*sparse( n_add .+ vely_ind, vort_ind,
        ones(size(vort_ind)), nrows, ncols);

    return W
end


function Lap_eigs( grid::T ) where T <: Grid
    # eigenvalues of RC (negative of the evals of the 5point stencil Lap)
    m = grid.nx;
    n = grid.ny;

    ii, jj = meshgrid( 1:(m-1), 1:(n-1) );  # Should replace with broadcasting
    Λ = -2*( cos.( π*ii/m ) .+ cos.( π*jj/n ) .- 2);
    return Λ
end


function Λinv_fn!(x::AbstractArray,
                  b::AbstractArray,
                  m::Int,
                  n::Int,
                  Λ::AbstractArray,
                  dst_plan::Tuple{Any, Array{Float64, 2}})
    """
    Solve inverse Laplacian RC (or similar for A matrix)
    Used to construct both Ainv and RCinv operators
    """
    # reshape for inversion in fourier space
    b = reshape( b, m-1, n-1)
    x = reshape( x, m-1, n-1)
    dst_inv!(x, b, Λ, dst_plan);
    # Include scale to make fwd/inv transforms equal
    rmul!(x, 4.0/( m*n ))
    x = reshape( x, (m-1)*(n-1), 1 )
end




function get_RCinv( grid::T,
                    Λ::AbstractArray,
                    dst_plan::Tuple{Any, Array{Float64, 2}}) where T <: Grid
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Benchmarks with nΓ = 39601
    Original function (non-mutating DST):
        1.506 ms (10 allocations: 928.73 KiB)
    LinearMap (mutating)
        1.251 ms (4 allocations: 160 bytes)
    """

    # give output in same size as input b (before being reshaped)
    return LinearMap((x, b) -> Λinv_fn!(x, b, grid.nx, grid.ny, Λ, dst_plan),
                     grid.nΓ; issymmetric=true, ismutating=true)
end
