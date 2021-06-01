#Driver function that steps forward in time
function advance!(t::Float64,
                  state::IBState{UniformGrid},
                  prob::IBProblem)
    grid = prob.model.grid

    # get irrotational flow contribution
    base_flux!(t, state, grid, prob.model.freestream)

    if MotionType(prob.model.bodies) != Static
        update_coupling!(prob.model, t)
        prob.Binv = get_B(prob.model, prob.Ainv)
    end

    # Alias working memory for notational clarity
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state


    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)
    circ2_st_vflx!( state.ψ, state.q, state.Γ, prob.model );

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    state.cfl = maximum( abs.( (1/grid.h^2) * state.q * dt ) ) ;
end


function advance!(t::Float64,
                  state::IBState{MultiGrid},
                  prob::IBProblem)
    """
    After updating get_nonlin:
        84.443 ms (1961 allocations: 76.93 MiB)
    After optimization:
        40.325 ms (719 allocations: 15.66 MiB)
    """
    grid = prob.model.grid

    # get irrotational flow contribution
    base_flux!(t, state, grid, prob.model.freestream)

    if MotionType(prob.model.bodies) != Static
        update_coupling!(prob.model, t)
        prob.Binv = get_B(prob.model, prob.Ainv)
    end

    # Alias working memory for notational clarity
    #   This leaves work.Γ2, Γ3 and work.q1, q2 available
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state


    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)
    circ2_st_vflx!( state.ψ, state.q, state.Γ, prob.model, grid.mg)

    # --Update circulation on all grids based on fine-grid correction
    for lev = 2:grid.mg
        @views coarsify!( state.Γ[:, lev-1], state.Γ[:, lev], grid)
    end

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    state.cfl = maximum( abs.( (1/(grid.h^2)) * @view(state.q[:, 1]) * dt ) )
end



#Back out background (freestream) flux q0 that is irrotational
function base_flux!(t::Float64,
                    state::IBState{UniformGrid},
                    grid::UniformGrid,
                    freestream::NamedTuple)
        """
        Initialize irrotational freestream flux
        """
    m = grid.nx;
    n = grid.ny;
    Ux = freestream.Ux(t)
    Uy = freestream.Uy(t)
    α = freestream.inclination(t)
    state.q0[ 1:(m-1)*n ] .= (Ux*cos(α) - Uy*sin(α))* grid.h  # x-flux
    state.q0[ (m-1)*n+1:end ] .= (Ux*sin(α) + Uy*cos(α))* grid.h  # y-flux
end


function base_flux!(t::Float64,
                    state::IBState{MultiGrid},
                    grid::MultiGrid,
                    freestream::NamedTuple)
        """
        Initialize irrotational freestream flux
        """
    m = grid.nx;
    n = grid.ny;
    Ux = freestream.Ux(t)
    Uy = freestream.Uy(t)

    α = freestream.inclination(t)
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:(m-1)*n, lev ] .= (Ux*cos(α) - Uy*sin(α))* hc    # x-flux
        state.q0[ (m-1)*n+1:end, lev ] .= (Ux*sin(α) + Uy*cos(α))*hc # y-flux
    end
end


function get_trial_state!(qs::AbstractArray,
                        Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    Combine explicit Laplacian and nonlinear terms into a rhs
       then use Ainv to return trial circulation Γs

    High-level version of AB2:
    rhs = A*Γ .-
    #    3*dt/2 * nonlin .+
    #    dt/2 * nonlin_prev

    Then do Ainv of that to back out trial circ
    """
    dt = prob.scheme.dt
    work = prob.work
    rhs = work.Γ2  # RHS of discretized equation

    #compute the nonlinear term for the current time step
    get_nonlin!( state.nonlin[1], state, prob );

    # Explicit part of Laplacian
    mul!(rhs, prob.A, state.Γ)

    # Explicit nonlinear terms from multistep scheme
    for n=1:length(prob.scheme.β)
        work.Γ3 .= state.nonlin[n]
        rmul!(work.Γ3, prob.scheme.β[n]*dt)
        rhs .-= work.Γ3
    end

    # Store current nonlinear term
    state.nonlin[2] .= state.nonlin[1];

    # Trial circulation  Γs = Ainv * rhs
    mul!(Γs, prob.Ainv, rhs);

    # Trial velocity  (note ψ is used here as a dummy variable)
    circ2_st_vflx!( state.ψ, qs, Γs, prob.model );
end



function get_trial_state!(qs::AbstractArray,
                        Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    Combine explicit Laplacian and nonlinear terms into a rhs
       then use Ainv to return trial circulation Γs

    High-level version of AB2:
    rhs = A*Γ .-
    #    3*dt/2 * nonlin .+
    #    dt/2 * nonlin_prev .+
    #    dt/2 * rhsbc

    Then do Ainv of that to back out trial circ
    """
    dt = prob.scheme.dt
    work = prob.work
    grid = prob.model.grid
    rhsbc = work.bc
    rhs = work.Γ2  # RHS of discretized equation

    for lev = grid.mg:-1:1

        # compute the nonlinear term for the current time step
        get_nonlin!( @view(state.nonlin[1][:, lev]), state, prob, lev );

        # contribution of Laplacian term...
        rhsbc .*= 0.0
        if lev < grid.mg
            # from explicit treatment of circulation
            get_Lap_BCs!( rhsbc, @view(state.Γ[:, lev+1]), lev, prob.model )

            # from current (trial) vorticity at previous grid level
            get_Lap_BCs!( rhsbc, @view(Γs[:,lev+1]), lev, prob.model)
        end

        # Combine explicit Laplacian and nonlinear terms into a rhs
        @views mul!(rhs[:, lev], prob.A[lev], state.Γ[:, lev])

        for n=1:length(prob.scheme.β)
            work.Γ3[:, lev] .= state.nonlin[n][:, lev]
            rmul!(work.Γ3, prob.scheme.β[n]*dt)
            rhs .-= work.Γ3
        end

        # Include boundary conditions
        #   High-level: rhs += 0.5*rhsbc
        work.Γ3[:, lev] .= rhsbc
        rmul!(work.Γ3, 0.5*dt)
        rhs .+= work.Γ3

        # Trial circulation  Γs = Ainv * rhs
        # TODO: use @view to do in-place multiplication
        # Doesn't work here because of view and FFT plan for even indices... WHY??
        Γs[:, lev] .= prob.Ainv[lev] * rhs[:, lev]
    end

    # Store nonlinear solution for use in next time step
    state.nonlin[2] .= state.nonlin[1]

    # Trial velocity on 1st grid level (don't need all grids)
    #   Streamfunction state.ψ is a dummy variable here to compute qs from Γs
    grid = prob.model.grid
    circ2_st_vflx!( state.ψ, qs, Γs, prob.model, 2 );
end


function boundary_forces!(F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    """
    Solve the modified Poisson equation (26)
    Dispatch based on the type of motion in the problem
    Allows precomputing regularization and interpolation for Static motions
    """
    boundary_forces!( MotionType(prob.model.bodies), F̃b, qs, q0, prob)
end

function boundary_forces!(::Type{Static},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    """
    Solve the Poisson equation (25) in Colonius & Taira (2008)
        for uB = 0 and bc2 = 0
        Bf̃ = Eq
           = ECψ
    """
    E = prob.model.mats.E
    h = prob.model.grid.h
    # Working memory for in-place operations
    qwork = @view(prob.work.q2[:, 1])
    broadcast!(+, qwork, qs, q0)                     # qs + q0
    mul!(F̃b, E, qwork)                               # E*(qs .+ state.q0)... using fb here as working array
    F̃b .= (1/h)*prob.Binv*F̃b                         # Allocates a small amount of memory
end



function boundary_forces!(::Type{RotatingCyl},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    """
    Solve the Poisson equation (25) in Colonius & Taira (2008)
        for bc2 = 0 with specialized situation of rotating cylinder
    In this case the points don't need to move, but do have nonzero velocity
        Bf̃ = Eq - ub
           = ECψ - ub
    """
    E = prob.model.mats.E
    h = prob.model.grid.h

    # Working memory for in-place operations
    F̃work = similar(F̃b)
    qwork = @view(prob.work.q2[:, 1])

    broadcast!(+, qwork, qs, q0)           # qs + q0
    mul!(F̃work, E, qwork)                  # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(F̃b, prob.Binv, F̃work/h);
    #F̃b .= (1/h)*prob.Binv*F̃b               # Allocates a small amount of memory
end



function project_circ!(Γs::AbstractArray,
                     state::IBState,
                     prob::IBProblem)
    """
    Update circulation to satisfy no-slip condition

    Dispatch based on the type of motion in the problem
    Allows precomputing regularization and interpolation for Static motions
    """
    project_circ!(MotionType(prob.model.bodies), Γs, state, prob)
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- Array(prob.Ainv[1] * (mats.RET*fb_til_dt) )
    """
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    RET = prob.model.mats.RET   # Precomputed R * E'
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs   # Now Γs is free for working memory
    @views mul!( Γs[:, 1], RET, fb_til_dt)
    @views mul!( Γwork, prob.Ainv, Γs[:, 1])
    state.Γ[:, 1] .-= Γwork
end


function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    High-level version:s
        state.Γ[:, 1] .= Γs .- Array(prob.Ainv[1] * (mats.RET*fb_til_dt) )
    """
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    RET = prob.model.mats.RET   # Precomputed R * E'
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], RET, fb_til_dt)
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])  # This is the only difference with the UniformGrid version
    state.Γ[:, 1] .-= Γwork
end

#Utilities for storing stress values
function update_stress!(state::IBState,
                        prob::IBProblem)
    """
    Store surface stresses and integrated forces

    Mutates "state"
    """
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    fb_til_dt = state.F̃b

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    state.CD = [state.CD; zeros(length(prob.model.bodies))]
    for j = 1 : length(prob.model.bodies)
        ds = prob.model.bodies[j].ds
        # surface stresses
        state.fb[j] .= fb_til_dt[nbod_tally .+ (1:nf[j])] *(h / dt) ./ [ds; ds]

        # integrated forces
        state.CD[size(state.CD,1),j] = 2 * sum( ds .* state.fb[j][1 : nb[j]] )
        state.CL[j] = 2 * sum( ds .* state.fb[j][1 .+ nb[j] : nf[j] ] )

        # update body index
        nbod_tally += nf[j];
    end
end


function AB2(dt::Float64)
    """
    Initialize second-order Adams-Bashforth scheme
    """
    return AdamsBashforth(dt, [1.5, -0.5])
end




"""
For dispatching to Static motions
function static_fn(model::IBModel{<:Grid, <:Body{Static}})
For other Motions
function dynamic_fn(model::IBModel)
"""

function update_coupling!(model::IBModel, t::Float64)
    bodies, grid = model.bodies, model.grid
    for j=1:length(bodies)
        move_body!(bodies[j], t)
    end

    model.mats.E = coupling_mat( grid, bodies )
    model.mats.RET = (model.mats.E*model.mats.C)'
end
