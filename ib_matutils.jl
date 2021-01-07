"""
Utility functions for constructing sparse IBPM matrices
"""

# Duplicate of Matlab function
repelem( x, N ) = Array( reshape( repeat(x, 1, N)', : ) );

function meshgrid(xin::Union{Array{Number, 1}, UnitRange{Int64}},
                  yin::Union{Array{Number, 1}, UnitRange{Int64}})
    """
    duplicate of matlab meshgrid function
    Should be able to get rid of this with broadcasting
    """
    nx=length(xin)
    ny=length(yin)
    xout=zeros(ny,nx)
    yout=zeros(ny,nx)
    for jx=1:nx
        for ix=1:ny
            xout[ix,jx]=xin[jx]
            yout[ix,jx]=yin[ix]
        end
    end
    return (x=xout, y=yout)
end

function get_velx_ind( i::Int, j::Int, grid::UniformGrid )
    """
    Return vector indices for x - velocity (velx_ind)
    corresponding to 2d array indices.

     e.g. get_velx_ind(2, 3, grid) returns 2*(m-1) + 2, which is the
     vector index of the x velocity at the point with xindex = 2, yindex = 3,
     and gridlevel = 1.

    Note: x velocity on grid 1 is of size (m-1)*n.
       For all other grids it is of size (m-1)*n/2 + m/2*n/2.
    """
    m = grid.nx;
    n = grid.ny;

    # First grid level
    return (m-1) .* (j-1) + i ;
end


function get_velx_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::UniformGrid )
    """
    Return vector indices for x - velocity (velx_ind)
    corresponding to 2d array indices.
    """
    @assert length(x_idx) == length(y_idx)
    return [get_velx_ind(x_idx[i], y_idx[i], grid) for i=1:length(x_idx)]
end


function get_velx_ind( i::Int, j::Int, grid::MultiGrid; glev::Int=1 )
    """
    Return vector indices for x - velocity (velx_ind)
    corresponding to 2d array indices.

     e.g. get_velx_ind(2, 3, grid) returns 2*(m-1) + 2, which is the
     vector index of the x velocity at the point with xindex = 2, yindex = 3,
     and gridlevel = 1.

    Note: x velocity on grid 1 is of size (m-1)*n.
       For all other grids it is of size (m-1)*n/2 + m/2*n/2.
    """
    m = grid.nx;
    n = grid.ny;

    # First grid level
    if glev == 1
        return (m-1) .* (j-1) + i
    else
        #  index starts from the end of the last gridlevel...
        n_add = (m-1) * n +  #  contribution from 1st gridlevel
            (glev - 2) * ( (m-1)*n/2 + m/2*n/2 ); # contribution from remaining glevs

        # points that are below finer grid:
        velx_ind = (j <= n/4) .* ( n_add + (m-1) .* (j-1) + i );

        # points that contain finer grid and are to left of finer grid:
        n_bott = (m-1)*n/4;
        velx_ind = velx_ind + ( j > n/4 & j <= 3*n/4 & i <= m/4 ) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + i );

        # points that contain finer grid and are to right of finer grid:
        velx_ind = velx_ind + (j > n/4 & j <= 3*n/4 & i >= 3*m/4) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + (i-m/2+1) );

        # points that are above finer grid:
        n_mid = m/2 * n/2;
        velx_ind = velx_ind + (j > 3*n/4 ).*
            ( n_add + n_bott + n_mid + (j-3*n/4-1)*(m-1) + i );

        return velx_ind
    end
end

function get_velx_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::MultiGrid; glev::Int=1 )
    """
    Return vector indices for x - velocity (velx_ind)
    corresponding to 2d array indices.

     e.g. get_velx_ind(2, 3, grid) returns 2*(m-1) + 2, which is the
     vector index of the x velocity at the point with xindex = 2, yindex = 3,
     and gridlevel = 1.

    Note: x velocity on grid 1 is of size (m-1)*n.
       For all other grids it is of size (m-1)*n/2 + m/2*n/2.
    """
    @assert length(x_idx) == length(y_idx)
    return [get_velx_ind(x_idx[i], y_idx[i], grid; glev=glev) for i=1:length(x_idx)]
end


function get_vely_ind( i::Int, j::Int, grid::UniformGrid )
    """
    %Return vector indices for y - velocity (vely_ind)
    %corresponding to 2d array indices.

    % e.g. get_vely_ind(2, 3, 1, grid_parms) returns 2*m + 2, which is the
    % vector index of the x velocity at the point with xindex = 2, yindex = 3,
    % and gridlevel = 1.

    %Note: x velocity on grid 1 is of size m*(n-1).
    %   For all other grids it is of size m*n/2 + m/2*(n/2-1).
    """
    m = grid.nx;
    n = grid.ny;

    # First grid level
    return m .* (j-1) + i ;
end


function get_vely_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::UniformGrid )
    @assert length(x_idx) == length(y_idx)
    return [get_vely_ind(x_idx[i], y_idx[i], grid) for i=1:length(x_idx)]
end


function get_vely_ind( i::Int, j::Int, grid::MultiGrid; glev::Int=1 )
    """
    %Return vector indices for y - velocity (vely_ind)
    %corresponding to 2d array indices.

    % e.g. get_vely_ind(2, 3, 1, grid_parms) returns 2*m + 2, which is the
    % vector index of the x velocity at the point with xindex = 2, yindex = 3,
    % and gridlevel = 1.

    %Note: x velocity on grid 1 is of size m*(n-1).
    %   For all other grids it is of size m*n/2 + m/2*(n/2-1).
    """
    m = grid.nx;
    n = grid.ny;

    # 1st gridlevel is easy...
    if glev == 1
        return m .* (j-1) + i

    # If not 1st grid level...
    else
        # index starts from the end of the last gridlevel...
        n_add = m * (n-1) +  # contribution from 1st gridlevel
            (glev - 2) * ( m*n/2 + m/2*(n/2-1) ); # contribution from remaining glevs


        # points that are below finer grid:
        vely_ind = (j <= n/4) .* ( n_add + m .* (j-1) + i );

        # points that contain finer grid and are to left of finer grid:
        n_bott = m*n/4;
        vely_ind = vely_ind + ( j > n/4 & j < 3*n/4 & i <= m/4 ) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + i );

        # points that contain finer grid and are to right of finer grid:
        vely_ind = vely_ind + (j > n/4 & j < 3*n/4 & i > 3*m/4) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + (i-m/2) );

        # points that are above finer grid:
        n_mid = m/2 * (n/2-1);
        vely_ind = vely_ind + (j >= 3*n/4 ).*
            ( n_add + n_bott + n_mid + (j-3*n/4)*m + i );
        return vely_ind

    end

end


function get_vely_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::MultiGrid; glev::Int=1 )
    @assert length(x_idx) == length(y_idx)
    return [get_vely_ind(x_idx[i], y_idx[i], grid; glev=glev) for i=1:length(x_idx)]
end



function get_vort_ind( i::Int, j::Int, grid::UniformGrid )
    """
    %Return vector indices for vorticity (vort_ind)
    %corresponding to 2d array indices.

    % e.g. get_vort_ind(2, 3, 1, grid_parms) returns 2*(m-1) + 2, which is the
    % vector index of the x velocity at the point with xindex = 2, yindex = 3,
    % and gridlevel = 1.

    %Note: vorticity on grid 1 is of size (m-1)*(n-1).
    %   For all other grids it is of size (m-1)*n/2 + m/2*(n/2-1).
    """
    m = grid.nx;
    n = grid.ny;

    # First grid level
    return (m-1) .* (j-1) + i ;
end


function get_vort_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::UniformGrid )
    @assert length(x_idx) == length(y_idx)
    return [get_vort_ind(x_idx[i], y_idx[i], grid) for i=1:length(x_idx)]
end


function get_vort_ind( i::Int, j::Int, grid::MultiGrid; glev::Int=1 )
    """
    %Return vector indices for vorticity (vort_ind)
    %corresponding to 2d array indices.

    % e.g. get_vort_ind(2, 3, 1, grid_parms) returns 2*(m-1) + 2, which is the
    % vector index of the x velocity at the point with xindex = 2, yindex = 3,
    % and gridlevel = 1.

    %Note: vorticity on grid 1 is of size (m-1)*(n-1).
    %   For all other grids it is of size (m-1)*n/2 + m/2*(n/2-1).
    """
    m = grid.nx;
    n = grid.ny;

    # 1st gridlevel is easy...
    if glev == 1
        return (m-1) .* (j-1) + i

    # If not 1st grid level...
    else
        # index starts from the end of the last gridlevel...
        n_add = (m-1) * (n-1) +      # contribution from 1st gridlevel
            (glev - 2) * ( (m-1)*n/2 + m/2*(n/2-1) ); # contribution from remaining glevs


        # points that are below finer grid:
        vort_ind = (j <= n/4) .* ( n_add + (m-1) .* (j-1) + i );

        # points that contain finer grid and are to left of finer grid:
        n_bott = (m-1)*n/4;
        vort_ind = vort_ind + (j > n/4 & j < 3*n/4 & i <= m/4 ) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + i );

        # points that contain finer grid and are to right of finer grid:
        vort_ind = vort_ind + (j > n/4 & j < 3*n/4 & i >= 3*m/4) .*
            ( n_add + n_bott + m/2*(j-n/4-1) + (i-m/2+1) );

        # points that are above finer grid:
        n_mid = m/2*( n/2 - 1 );
        vort_ind = vort_ind + (j >= 3*n/4 ).*
            ( n_add + n_bott + n_mid + (j-3*n/4)*(m-1) + i );
        return vort_ind

    end

end

function get_vort_ind( x_idx::Array{Int, 1}, y_idx::Array{Int, 1}, grid::MultiGrid; glev::Int=1 )
    @assert length(x_idx) == length(y_idx)
    return [get_vort_ind(x_idx[i], y_idx[i], grid; glev=glev) for i=1:length(x_idx)]
end


function get_vort_bounds(grid::MultiGrid; glev::Int=1 )
    """
    Get indices on grid level (glev) that corresponds to the left,
    right, top, and bottom indices on the finer grid level.
    """
    m = grid.nx;
    n = grid.ny;

    # [left, right, bottom, top]
    return [get_vort_ind( m/4, n/4 : 3*n/4, grid; glev=glev )
            get_vort_ind( 3*m/4, n/4 : 3*n/4, grid; glev=glev )
            get_vort_ind( m/4 : 3*m/4, n/4, grid; glev=glev )
            get_vort_ind( m/4 : 3*m/4, 3*n/4, grid; glev=glev )]

end
