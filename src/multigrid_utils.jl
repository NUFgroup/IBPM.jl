

function coarsify!( circf::AbstractVector,
                    circc::AbstractVector,
                    grid::MultiGrid;
                    verbose::Bool=false )
    """
    %Take circulation on a finer grid (circf) and use it to replace the
    %circulation on a coarser grid (circc) in the overlapping region.

    Original:
        2.417 ms (167 allocations: 4.05 MiB)
    Optimized:
        173.699 μs (8 allocations: 192 bytes)
    """
    m = grid.nx;
    n = grid.ny;

    addc = (n÷4)*(m-1) + (m÷4)

    for i=0:(m÷2-1)*(n÷2-1)-1
        addf = 2*(m-1)*(i÷(m÷2-1))
        ridx = mod(i, m÷2-1)  # "Row-like" index

        # Overlap for coarse grid
        coarse_c = addc + 1+ridx + addf÷2

        # Fine grid points
        fine_c = addf + (m-1) + 2*(1+ridx)  # Center point

        fine_t = addf + 2*(m-1) + 2*(1+ridx)       # Top point
        fine_b = addf + 2*(1+ridx)                 # Bottom point
        fine_l = addf + (m-1) + 2*(1+ridx)-1       # Left point
        fine_r = addf + (m-1) + 2*(1+ridx)+1       # Right point
        fine_tl = addf + 2*(m-1) + 2*(1+ridx)-1    # Top left point
        fine_tr = addf + 2*(m-1) + 2*(1+ridx)+1    # Top right point
        fine_bl = addf + 2*(1+ridx)-1              # Bottom left point
        fine_br = addf + 2*(1+ridx)+1              # Bottom right point

        #--interpolate finer circulation onto coarser one
        circc[coarse_c] = circf[ fine_c ] +
            0.5 * ( circf[fine_t] + circf[fine_b] +
                    circf[fine_l] + circf[fine_r] ) +
            0.25 * ( circf[fine_tl] + circf[fine_tr] +
                     circf[fine_bl] + circf[fine_br] )
    end



end

function get_stfn_BCs!( Γb::AbstractArray,
                        Γ::AbstractArray,
                        model::IBModel{MultiGrid, <:Body} )
    """
    BCs used for RCinv (taken from the streamfunction on the coarser grid) and curl

    Not used for single-grid

        # bc: [bottom, top, left, right]
    Note: same indexing as in get_Lap_BCs
    """
    grid = model.grid
    m = grid.nx;
    n = grid.ny;
    mg = grid.mg;
    m = grid.nx
    n = grid.ny

# Bottom part
    naddc = (n÷4 - 1)*(m-1) + m÷4  # Add to coarse indices
    naddf = 1   # Add to fine indices

    # points that need to average coarser domain:
    for i=0:m÷2-1
        Γb[naddf+2*i, 1] += 0.5*( Γ[naddc+i] + Γ[naddc+i+1] )
    end

    # points that don't need to average coarser domain:
    for i=1:m÷2-1
        Γb[naddf+(2*i-1), 1] += Γ[naddc+i]
    end

# Top part
    naddc = (3*n÷4 - 1)*(m-1) + m÷4  # Add to coarse indices
    naddf = (n-2)*(m-1) + 1  # Add to fine indices

    # indices on coarse grid corresponding to top part of fine grid
    for i=0:m÷2-1
        Γb[naddf + 2*i, 2] += 0.5*( Γ[naddc+i] + Γ[naddc+i+1] )
    end

    # points that don't need to average coarser domain:
    for i=1:m÷2-1
        Γb[naddf+(2*i-1), 2] += Γ[naddc+i]
    end

# left part
    naddc = (n÷4 - 1)*(m-1) + m÷4
    naddf = 1;

    # indices on coarse grid corresponding to left edge of fine grid
    for i=0:n÷2-1
        Γb[naddf + 2*i*(m-1), 3] += 0.5*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
    end

    # points that don't need to average coarser domain:
    for i=1:n÷2-1
        Γb[naddf+(2*i-1)*(m-1), 3] += Γ[naddc+i*(m-1)]
    end

# right part
    naddc = (n÷4 - 1)*(m-1) + 3*m÷4
    naddf = m-1;

    # indices on coarse grid corresponding to right edge of fine grid
    for i=0:n÷2-1
        Γb[naddf + 2*i*(m-1), 4] += 0.5*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
    end

    # points that don't need to average coarser domain:
    for i=1:n÷2-1
        Γb[naddf+(2*i-1)*(m-1), 4] += Γ[naddc+i*(m-1)]
    end
# **

end


function curl!( q::AbstractArray,
                ψ::AbstractArray,
                stbc::AbstractArray,
                model::IBModel{MultiGrid, <:Body} )
    """
    %Compute velocity flux from streamfunction.
    %   Note: requires streamfunction from coarser grid on edge
    %         of current domain (stored in stbc)

    stbc = [bottom, top, left, right]
    """
    m = model.grid.nx;
    n = model.grid.ny;

    # First get contribution without accounting for BCs:
    mul!(q, model.mats.C, ψ)  # q = C*ψ

    #!!x-velocity block (terms on bottom and top edges):
        ### Bottom part
        q[ 1:m-1 ] .-= stbc[ 1:m-1, 1 ] ;

        ### Top part
        naddc = (n-1)*(m-1)
        naddf = (n-2)*(m-1)
        for i=1:m-1
            q[naddc + i] += stbc[naddf + i, 2]
        end

    #!!y-velocity block
        ### Left part
        nadd = (m-1)*n;  # yvel index begins after x vels
        for i=0:n-2
            q[nadd + 1+m*i] += stbc[1+i*(m-1), 3]
        end

        ### right part
        for i=1:n-1
            q[nadd + m*i] -= stbc[i*(m-1), 4]
        end
end


function get_Lap_BCs!( rhsbc::AbstractArray, rhs::AbstractArray, lev::Int,
                        model::IBModel{MultiGrid, <:Body} )
    """
    BCs used for Ainv in multigrid

    Original:
        38.069 μs (64 allocations: 352.59 KiB)
    Optimized:
        27.111 μs (3 allocations: 309.53 KiB)
    """

    # scaling factor to multiply gamma by:
    #   The 1/4 turns circ on coarser grid into circ on finer grid
    #   The 1/hc^2 is because the Laplacian is a second deriv operator
    #   The 1/Re number is the physical scaling for the viscous Lap
    hc = model.grid.h * 2^(lev-1);  # grid spacing on current grid
    scl = ( 1/4 / (hc^2) ) / model.Re;

    # Call function to loop over boundary points and average from coarser grid
    rhs_BCs!(rhsbc, rhs, model.grid, scl)
end



function rhs_BCs!( Γb::AbstractArray, Γ::AbstractArray,
                  grid::MultiGrid, scl::Float64)
    """
    Function to abstract common indexing loops in multi-grid boundary conditions
    Γb:  (nΓ x 1) fine-grid vector that will get edited
    Γ:   (nΓ x 1) vector on coarser grid

    TO DO: CAN THIS BE MERGED WITH WΓ_BCs and Qqq0_BCs in get_nonlinear and stfn_BCs above??
    """
    m = grid.nx
    n = grid.ny

    # Bottom part
        naddc = (n÷4 - 1)*(m-1) + m÷4  # Add to coarse indices
        naddf = 1   # Add to fine indices

        # points that need to average coarser domain:
        for i=0:m÷2-1
            Γb[naddf+2*i] += 0.5*scl*( Γ[naddc+i] + Γ[naddc+i+1] )
        end

        # points that don't need to average coarser domain:
        for i=1:m÷2-1
            Γb[naddf+(2*i-1)] += scl*Γ[naddc+i]
        end

    # Top part
        naddc = (3*n÷4 - 1)*(m-1) + m÷4  # Add to coarse indices
        naddf = (n-2)*(m-1) + 1  # Add to fine indices

        # indices on coarse grid corresponding to top part of fine grid
        for i=0:m÷2-1
            Γb[naddf + 2*i] += 0.5*scl*( Γ[naddc+i] + Γ[naddc+i+1] )
        end

        # points that don't need to average coarser domain:
        for i=1:m÷2-1
            Γb[naddf+(2*i-1)] += scl*Γ[naddc+i]
        end

    # left part
        naddc = (n÷4 - 1)*(m-1) + m÷4
        naddf = 1;

        # indices on coarse grid corresponding to left edge of fine grid
        for i=0:n÷2-1
            Γb[naddf + 2*i*(m-1)] += 0.5*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
        end

        # points that don't need to average coarser domain:
        for i=1:n÷2-1
            Γb[naddf+(2*i-1)*(m-1)] += scl*Γ[naddc+i*(m-1)]
        end

    # right part
        naddc = (n÷4 - 1)*(m-1) + 3*m÷4
        naddf = m-1;

        # indices on coarse grid corresponding to right edge of fine grid
        for i=0:n÷2-1
            Γb[naddf + 2*i*(m-1)] += 0.5*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
        end

        # points that don't need to average coarser domain:
        for i=1:n÷2-1
            Γb[naddf+(2*i-1)*(m-1)] += scl*Γ[naddc+i*(m-1)]
        end
    # **

end
