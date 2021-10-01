using LinearAlgebra: norm  # FOR DEBUGGING

"""
    advance!(state::IBState, prob::IBProblem, t::Float64)

    Advance state forward in time.
"""
function advance!(state::IBState{MultiGrid},
                  prob::IBProblem,
                  t::Float64)
    grid = prob.model.grid
    #Update either the grid motion or body motion based on user-prescribed
    #kinematics. There can be additional FSI on top of this.
    update_prescribed_body_motion!(state, prob, t)

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
    # @views update_surface_quantities!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)

    @show maximum(abs.(state.fb)), maximum(abs.(state.F̃b))
    @views update_surface_quantities!(qs[:, 1], state, prob)
    @show maximum(abs.(state.fb)), maximum(abs.(state.F̃b))
    update_stress!(state, prob) #Compute integral quantities and store in state
    @show maximum(abs.(state.fb)), maximum(abs.(state.F̃b))

    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)

    # Interpolate values from finer grid to center region of coarse grid
    vort2flux!( state.ψ, state.q, state.Γ, prob.model, grid.mg )

    return nothing
end
