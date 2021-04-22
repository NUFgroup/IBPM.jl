using LinearAlgebra: norm  # FOR DEBUGGING

"""
    advance!(state::IBState, prob::IBProblem, t::Float64)

Advance state forward in time.
"""
function advance!(state, prob, t) end

function advance!(state::IBState{UniformGrid},
                  prob::IBProblem,
                  t::Float64)
    grid = prob.model.grid

    # Move bodies and update coupling matrices (E)
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == MovingGrid
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    qs = prob.model.work.q1  # Trial flux
    Γs = prob.model.work.Γ1  # Trial circulation

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
    vort2flux!( state.ψ, state.q, state.Γ, prob.model );

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    state.cfl = maximum( abs.( (1/grid.h^2) * state.q * prob.scheme.dt ) ) ;

    return nothing
end



function advance!(state::IBState{MultiGrid},
                  prob::IBProblem,
                  t::Float64)
    grid = prob.model.grid
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == MovingGrid
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    #   This leaves work.Γ2, Γ3 and work.q1, q2 available
    qs = prob.model.work.q1  # Trial flux
    Γs = prob.model.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state

    #println(sum(Γs.^2))
    #println(sum(state.Γ[:, 1].^2))
    #println(sum(qs.^2))

    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)
    #println("Final circulation: ", sum(state.Γ.^2))

    # Interpolate values from finer grid to center region of coarse grid
    vort2flux!( state.ψ, state.q, state.Γ, prob.model, grid.mg );

    #println("Final circulation: ", sum(state.Γ.^2))
    #println("Final flux: ", sum(state.q.^2))
    #println("Final stfn: ", sum(state.ψ.^2))

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    # TODO: DOES THIS ALLOCATE??
    state.cfl = maximum( @. abs( (1/(grid.h^2))*state.q[:, 1]*dt ) ) ;

    return nothing
end


"""
    get_trial_state!(qs, Γs, state, prob)

Compute trial circulation Γs that doesn't satisfy no-slip BCs

Combine explicit Laplacian and nonlinear terms into a rhs
   then invert implicit part to return trial circulation Γs

High-level version of AB2:
rhs = A*Γ .-
      3*dt/2 * nonlin .+
      dt/2 * nonlin_prev .+
      dt/2 * rhsbc

Then do Ainv of that to back out trial circ
"""
function get_trial_state!(qs, Γs, state, prob) end

function get_trial_state!(qs::AbstractArray,
                          Γs::AbstractArray,
                          state::IBState{UniformGrid},
                          prob::IBProblem)
    grid = prob.model.grid  # SHOULDN'T NEED THIS HERE
    dt = prob.scheme.dt
    work = prob.model.work
    rhs = work.Γ2  # RHS of discretized equation

    #compute the nonlinear term for the current time step
    nonlinear!( state.nonlin[1], state, prob );

    # Explicit part of Laplacian
    @views mul!(rhs, prob.A, state.Γ[:, 1])

    # Explicit nonlinear terms from multistep scheme
    for n=1:length(prob.scheme.β)
        rhs .+= prob.scheme.β[n]*dt*@view(state.nonlin[n][:, 1]);
    end

    # Trial circulation  Γs = Ainv * rhs
    mul!(Γs, prob.Ainv, rhs);

    # Trial velocity  (note ψ is used here as a dummy variable)
    vort2flux!( state.ψ, qs, Γs, prob.model )

    # Store current nonlinear term
    state.nonlin[2] .= state.nonlin[1];

    return nothing
end

function get_trial_state!(qs::AbstractArray,
                          Γs::AbstractArray,
                          state::IBState{MultiGrid},
                          prob::IBProblem)
    dt = prob.scheme.dt
    grid = prob.model.grid
    rhsbc = prob.model.work.rhsbc
    #rhs = @view(work.Γ2[:, 1])  # RHS of discretized equation
    rhs = prob.model.work.Γ2  # RHS of discretized equation
    bc = prob.model.work.Γbc

    for lev=grid.mg:-1:1
        bc .*= 0.0; rhsbc .*= 0.0
        hc = grid.h * 2^( lev - 1);

        if lev < grid.mg
            @views get_bc!(bc, state.Γ[:, lev+1], grid)

            fac = 0.25*dt/ ( prob.model.Re * hc^2 )
            apply_bc!( rhsbc, bc, fac, grid )
        end

        #compute the nonlinear term for the current time step
        @views nonlinear!( state.nonlin[1][:, lev], state, bc, lev, prob );

        @views mul!( rhs, prob.A[lev], state.Γ[:, lev] )

        for n=1:length(prob.scheme.β)
            rhs .+= (prob.scheme.β[n]*dt)*@view(state.nonlin[n][:, lev])
        end

        # Include boundary conditions
        rhs .+= rhsbc

        # Trial circulation  Γs = Ainv * rhs
        @views mul!(Γs[:, lev], prob.Ainv[lev], rhs)
    end

    # Store nonlinear solution for use in next time step
    state.nonlin[2] .= state.nonlin[1]

    vort2flux!( state.ψ, qs, Γs, prob.model, grid.mg )
    return nothing
end

"""
    boundary_forces!(F̃b, qs, q0, prob)

Solve the Poisson equation (25) in Colonius & Taira (2008).

Dispatch based on the type of motion in the problem - allows precomputing
    regularization and interpolation where possible.
"""
function boundary_forces!(F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    boundary_forces!( MotionType(prob.model.bodies), F̃b, qs, q0, prob)
end

"""
    boundary_forces!(::Union{Type{Static}, Type{MovingGrid}},
                     F̃b, qs, q0, prob)

Solve modified Poisson problem for uB = 0 and bc2 = 0
```
 Bf̃ = Eq = ECψ
```
"""
function boundary_forces!(::Union{Type{Static}, Type{MovingGrid}},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    E = prob.model.mats.E
    h = prob.model.grid.h
<<<<<<< HEAD

    #qwork = @view(prob.model.work.q2[:, 1])  # Working memory for in-place operations
    Q = prob.model.work.q2  # Net flux

    broadcast!(+, Q, qs, q0)                     # qs + q0
    mul!(F̃b, E, Q)                               # E*(qs .+ state.q0)... using fb here as working array
    F̃b .= prob.Binv*F̃b                         # Allocates a small amount of memory

    return nothing
=======
    # Working memory for in-place operations
    qwork = @view(prob.work.q2[:, 1])
    broadcast!(+, qwork, qs, q0)                     # qs + q0
    mul!(F̃b, E, qwork)                               # E*(qs .+ state.q0)... using fb here as working array
    F̃b .= (1/h)*prob.Binv*F̃b                         # Allocates a small amount of memory
    #F̃b .= (1/h)*(prob.Binv\F̃b)                      # USE WITH CHOLESKY
>>>>>>> main
end

"""
    boundary_forces!(::Type{RotatingCyl}, F̃b, qs, q0, prob)

Solve the Poisson problem for bc2 = 0 with special case of rotating cylinder.

In this case the points don't need to move, but they do have nonzero velocity
```
Bf̃ = Eq - ub
   = ECψ - ub
```
"""
function boundary_forces!(::Type{RotatingCyl},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    E = prob.model.mats.E
    h = prob.model.grid.h

    # Working memory for in-place operations (small allocation)
    F̃work = similar(F̃b)
    #qwork = @view(prob.model.work.q2[:, 1])
    Q = prob.model.work.q2   # Net flux

    broadcast!(+, Q, qs, q0)           # qs + q0
    mul!(F̃work, E, Q)                  # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(F̃b, prob.Binv, F̃work);

    return nothing
end

"""
    project_circ!(Γs, state, prob)

Update circulation to satisfy no-slip condition.

Dispatch based on the type of motion in the problem.

This allows precomputing regularization and interpolation where possible.
"""
function project_circ!(Γs::AbstractArray,
                     state::IBState,
                     prob::IBProblem)
    project_circ!(MotionType(prob.model.bodies), Γs, state, prob)
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- prob.Ainv[1] * (mats.RET*fb_til_dt)
    """
    #Γwork = @view(prob.model.work.Γ3[:, 1]) # Working memory
    Γwork = prob.model.work.Γ2 # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs   # Now Γs is free for working memory
    @views mul!( Γs[:, 1], (E*C)', fb_til_dt)  # Γ = ∇ x (E'*fb)
    @views mul!( Γwork, prob.Ainv, Γs[:, 1])
    state.Γ[:, 1] .-= Γwork

    return nothing
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- prob.Ainv[1] * (mats.RET*fb_til_dt)
    """
    #println("=== PROJECT CIRC ===")
    #Γwork = @view(prob.model.work.Γ3[:, 1]) # Working memory
    Γwork = prob.model.work.Γ2 # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], (E*C)', fb_til_dt)  # Γ = ∇ x (E'*fb)
    #println(sum(Γs[:, 1].^2))
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])  # This is the only difference with the UniformGrid version
    #println(sum(Γwork.^2))

    state.Γ[:, 1] .-= Γwork

    #println(sum(state.Γ.^2))

    return nothing
end

"""
    update_stress!(state, prob)

Store surface stresses and integrated forces.

Mutates "state"
"""
function update_stress!(state::IBState,
                        prob::IBProblem)
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    fb_til_dt = state.F̃b

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j = 1 : length(prob.model.bodies)
        ds = prob.model.bodies[j].ds
        # surface stresses
        state.fb[j] .= fb_til_dt[nbod_tally .+ (1:nf[j])] *(h / dt) ./ [ds; ds] ;

        # integrated forces
        state.CD[j] = 2 * sum( ds .* state.fb[j][1 : nb[j]] ) ;
        state.CL[j] = 2 * sum( ds .* state.fb[j][1 .+ nb[j] : nf[j] ] ) ;

        # update body index
        nbod_tally += nf[j];
    end

    return nothing
end

"""
    AB2(dt::Float64)

Initialize second-order Adams-Bashforth scheme.
"""
function AB2(dt::Float64)
    return AdamsBashforth(dt, [1.5, -0.5])
end


"""
    update_bodies!(prob, t)

Update the immersed bodies and coupling matrices (if applicable).

TODO: break out by multiple dispatch... but don't duplicate code
"""
function update_bodies!(prob::IBProblem, t::Float64)
    model = prob.model
    bodies, grid = prob.model.bodies, prob.model.grid
    motion = MotionType(bodies)

    if motion != Static
        for j=1:length(bodies)
            move_body!(bodies[j], t)
        end
    end

    # For arbitrary motion in an inertial frame, have to update operators
    if motion == MotionFunction
        model.mats.E = coupling_mat( grid, bodies )
        model.mats.RET = (model.mats.E*model.mats.C)'
        prob.Binv = get_B(model, prob.Ainv)
    end
end
