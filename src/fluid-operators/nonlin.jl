

#Does this version get used anywhere?
function get_nonlin( state::IBState{UniformGrid}, prob::IBProblem )
    """
    Build the nonlinear term. Without accounting for BCs, the answer is
    nonlin = mats.R * ( ( mats.W * gamma ) .* ( mats.Q * (q + q0) ) );
    """

    # --parameters used in this function
    grid = prob.model.grid
    mats = prob.model.mats
    m = grid.nx;
    n = grid.ny;

    #---Build Wgamma and Q(q + q0) without accounting for BCs

    # the 1/hc^2 term is to convert circ to vort
    Wgam = 1/(grid.h^2) * mats.W * sol.Γ;
    return mats.R * ( Wgam .*  (mats.Q* (state.q + state.q0)));
end

function get_nonlin!( nonlin::AbstractArray,
                      state::IBState{UniformGrid},
                      prob::IBProblem )
    """
    Build the nonlinear term. Without accounting for BCs, the answer is
    nonlin = mats.R * ( ( mats.W * gamma ) .* ( mats.Q * (q + q0) ) );

    work.Γ4 and work.q2-3 are free for use here
    """

    # --parameters used in this function
    grid = prob.model.grid
    mats = prob.model.mats

    # --alias working memory
    WΓ = prob.work.q2
    qq0 = prob.work.q3
    Qqq0 = prob.work.q4

    # --in-place operations to compute nonlinear term
    #Wgam = 1/(hc^2) * mats.W * sol.Γ;
    mul!(WΓ, mats.W, state.Γ);
    rmul!(WΓ, 1/grid.h^2);

    broadcast!(+, qq0, state.q, state.q0);  # qq0 = state.q .+ state.q0
    mul!(Qqq0, mats.Q, qq0)         # Qqq0 = Q * qq0;

    mul!( nonlin, mats.R, WΓ.*Qqq0 )
end


function get_nonlin!( nonlin::AbstractVector,
                      state::IBState{MultiGrid},
                      prob::IBProblem,
                      lev::Int )
    """

    Build the nonlinear term. Without accounting for BCs, the answer is
    nonlin = mats.R * ( ( mats.W * gamma ) .* ( mats.Q * (q + q0) ) );

    However, we need to modify the W*gamma term and the Q*(q + q0) term to
    account for circulation and velocity flux terms from the coarse grid on
    the fine grid. (Note that R does not need to be modified).

    work.Γ3 and work.q2-3 are free for use here

    Original:
        1.946 ms (456 allocations: 8.88 MiB)
    Optimized:
        1.389 ms (24 allocations: 1.52 MiB)
    """

    # --parameters used in this function
    grid = prob.model.grid
    mats = prob.model.mats

    # --alias working memory
    WΓ = @view(prob.work.q2[:, 1])
    Qqq0 = @view(prob.work.q2[:, 2])
    qq0 = prob.work.q3

    # Coarse grid spacing
    hc = grid.h * 2^( lev - 1);

    #---Build qs=W*Γ and Q(q + q0) without accounting for BCs
    # the 1/hc^2 term is to convert circ to vort

    #Wgam = 1/(hc^2) * mats.W * sol.Γ[:, lev];
    mul!(WΓ, mats.W, @view(state.Γ[:, lev]));
    rmul!(WΓ, 1/hc^2);

    @views broadcast!(+, qq0[:, lev], state.q[:, lev], state.q0[:, lev]);  # qq0 = state.q .+ state.q0
    mul!(Qqq0, mats.Q, @view(qq0[:, lev]))         # Qqq0 = Q * qq0;

    # TODO: is this the same as one of the other functions....?
    #--Now take BCs into account:
    if (lev < grid.mg)
        # ** Wgamma term
        #scaling factor to multiply gamma by:
        #   The 1/4 turns circ on coarser grid into circ on finer grid
        #   The 1/2 is because this term is part of an average
        #   The hc^2 term converts circ to vort
        scl = 1/2 * 1/4 / hc^2;
        WΓ_BCs!( WΓ, @view(state.Γ[:, lev+1]), grid, scl )


        #** Qqq0 term
        #scaling factor to multiply qq0 by:
        #   The 1/2 converts the velocity flux on the coarser grid to a
        #   flux on the finer grid
        scl = 1/2 ;
        @views broadcast!(+, qq0[:, lev+1], state.q[:, lev+1], state.q0[:, lev+1]);  # qq0 = state.q .+ state.q0
        Qqq0_BCs!( Qqq0, @view(qq0[:, lev+1]), grid, scl )

    end

    #nonlin .= mats.R * ( WΓ .*  Qqq0);
    WΓ .*= Qqq0
    mul!( nonlin, mats.R, WΓ)
end


function WΓ_BCs!(  Γb::AbstractArray, Γ::AbstractArray,
                       grid::MultiGrid, scl::Float64)
    """
    Function to abstract indexing loops in W*Γ boundary conditions
    Γb:  (nΓ x 1) fine-grid vector that will get edited
    Γ:   (nΓ x 1) vector on coarser grid

    TO DO: CAN THIS BE MERGED WITH rhs_BCs and Qqq0_BCs??
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
    naddf = (n-1)*(m-1)  + 1 # Add to fine indices

    # indices on coarse grid corresponding to top part of fine grid
    for i=0:m÷2-1
        Γb[naddf + 2*i] += 0.5*scl*( Γ[naddc+i] + Γ[naddc+i+1] )
    end

    # points that don't need to average coarser domain:
    for i=1:m÷2-1
        Γb[naddf+(2*i-1)] += scl*Γ[naddc+i]
    end

    #!!y-velocity block (contributions from left and right edges)
    # indices for y-vel start after x velocities
    nadd = (m-1)*n;


#!!y-velocity block (contributions from left and right edges)
# indices for y-vel start after x velocities
nadd = (m-1)*n;

# left part
    naddc = (n÷4 - 1)*(m-1) + m÷4
    naddf = nadd + 1;

    # indices on coarse grid corresponding to left edge of fine grid
    for i=0:n÷2-1
        Γb[naddf + 2*i*m] += 0.5*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
    end

    # points that don't need to average coarser domain:
    for i=1:n÷2-1
        Γb[naddf+(2*i-1)*m] += scl*Γ[naddc+i*(m-1)]
    end

# right part
    naddc = (n÷4 - 1)*(m-1) + 3*m÷4
    naddf = nadd+m;

    # indices on coarse grid corresponding to right edge of fine grid
    for i=0:n÷2-1
        Γb[naddf + 2*i*m] += 0.5*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
    end

    # points that don't need to average coarser domain:
    for i=1:n÷2-1
        Γb[naddf+(2*i-1)*m] += scl*Γ[naddc+i*(m-1)]
    end
end





function Qqq0_BCs!( Γb::AbstractArray, Γ::AbstractArray,
                    grid::MultiGrid, scl::Float64)
    """
    Function to abstract indexing loops in Q*(q + q0) boundary conditions
    Γb:  (nΓ x 1) fine-grid vector that will get edited
    Γ:   (nΓ x 1) vector on coarser grid

    TO DO: CAN THIS BE MERGED WITH rhs_BCs and WΓ_BCs??
    """
    m = grid.nx
    n = grid.ny

    nadd = (m-1)*n;

    # Bottom part
        naddc = nadd + (n÷4 - 1)*m + m÷4 + 1  # Add to coarse indices
        naddf = 0   # Add to fine indices

        # points that need to average coarser domain:
        for i=1:m÷2-1
            Γb[naddf+2*i] -= 0.25*scl*( Γ[naddc+i-1] + Γ[naddc+i] )
        end

        # points that don't need to average coarser domain:
        for i=0:m÷2-1
            Γb[naddf+(2*i+1)] -= 0.5*scl*Γ[naddc+i]
        end

    # Top part
        naddc = nadd + (3*n÷4 - 1)*m + m÷4 + 1  # Add to coarse indices
        naddf = (n-1)*(m-1)  # Add to fine indices

        # indices on coarse grid corresponding to top part of fine grid
        for i=1:m÷2-1
            Γb[naddf+2*i] -= 0.25*scl*( Γ[naddc+i-1] + Γ[naddc+i] )
        end

        # points that don't need to average coarser domain:
        for i=0:m÷2-1
            Γb[naddf+(2*i+1)] -= 0.5*scl*Γ[naddc+i]
        end

    # left part
        naddc = (n÷4-1)*(m-1) + m÷4
        naddf = nadd+1;

        # indices on coarse grid corresponding to left edge of fine grid
        for i=1:n÷2-1
            Γb[naddf + (2*i-1)*m] += 0.25*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
        end

        # points that don't need to average coarser domain:
        for i=0:n÷2-1
            Γb[naddf+(2*i)*m] += 0.5*scl*Γ[naddc+(i+1)*(m-1)]
        end

    # right part
        naddc = (n÷4-1)*(m-1) + 3*m÷4
        naddf = nadd + m;

        # indices on coarse grid corresponding to left edge of fine grid
        for i=1:n÷2-1
            Γb[naddf + (2*i-1)*m] += 0.25*scl*( Γ[naddc+i*(m-1)] + Γ[naddc+(i+1)*(m-1)] )
        end

        # points that don't need to average coarser domain:
        for i=0:n÷2-1
            Γb[naddf+(2*i)*m] += 0.5*scl*Γ[naddc+(i+1)*(m-1)]
        end
    # **

end


function get_Q( grid::T ) where T <: Grid
    """
    %Build averaging operator Q = [Qx; Qy]. The x-velocity block (Qx) takes y
    %velocities and averages them onto the x-velocity edges, and the y-velocity
    %block (Qy) takes x velocities and averages them onto y-velocity edges.

    %Note: The Q used in the timestepping protocol postmultiplies
    %Minv to the Q obtained here.

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

    # x-vels to top left of current y-vel point
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

    %Note that W is not the full W used to time advance. The W for time stepping
    %is M * the W computed here.

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