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
    return nothing
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    High-level version:
        Γ = Γs - Ainv * (E*C)'*F̃b
    """
    Γwork = prob.model.work.Γ2 # Working memory
    E, C = prob.model.mats.E, prob.model.mats.C

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], (E*C)', state.F̃b[:, 1])  # Γ = ∇ x (E'*fb)
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])

    state.Γ[:, 1] .-= Γwork

    return nothing
end
