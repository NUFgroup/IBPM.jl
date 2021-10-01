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
