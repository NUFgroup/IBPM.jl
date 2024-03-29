"""
Selective frequency damping
Implemented as a wrapper around the standard time-stepper

NOTE: THIS WAS WRITTEN FOR THE MATLAB-BASED VERSION
    AND HASN'T BEEN TESTED FOR THE NEW FORTRAN-BASED ONE
"""

# Only need circulation and a time history of rhs term
#    for the auxiliary filtered state
mutable struct SFDState <: State
        Γ::Array{Float64, 2}     # Circulation
        nonlin::Array{Array{Float64, 2}, 1}  # Memory of nonlinear terms
end

mutable struct SFDProblem <: AbstractIBProblem
    ib_prob::IBProblem
    model::IBModel
    scheme::ExplicitScheme
    A
    Ainv
    Binv
    state::SFDState  # Filtered state
    Δ::Float64    # Filter width
    χ::Float64    # Damping
end


function init_sfd(prob::IBProblem,
                  x::State,
                  Δ::Float64,
                  χ::Float64)
    """
    Initialize selective frequency damping
    Strategy is to alias the IBProblem fields so that methods
    can be called normally, but keep the ib_prob to avoid
    duplicating code
    """
    x̄ = SFDState(deepcopy(x.Γ), deepcopy(x.nonlin))
    return SFDProblem(prob,
                      prob.model,
                      prob.scheme,
                      prob.A,
                      prob.Ainv,
                      prob.Binv,
                      x̄, Δ, χ)
end

function advance!(state::IBState,
                  sfd::SFDProblem,
                  t::Float64)
    """
    Wrapper around the main time-stepper
    Note that since this updates the circulation after the KKT solver,
       it assumes that the body is fixed and also that the damped state
       satisfies the boundary conditions (which it should, since it's
       just a linear combination of the primary states).
    """
    β, dt = sfd.scheme.β, sfd.scheme.dt

    # Advance primary state
    advance!(state, sfd.ib_prob, t)

    # --- Update states ---

    # Compute filtered term χ*(Γ - Γ̄)
    broadcast!(-, sfd.state.nonlin[1], state.Γ, sfd.state.Γ)

    # Explicit time stepping
    for n=1:length(β)
        @. state.Γ -= (β[n]*dt*sfd.χ)*sfd.state.nonlin[n]
        @. sfd.state.Γ += (β[n]*dt/sfd.Δ)*sfd.state.nonlin[n]
    end

    # Update primary state streamfunction and velocity flux
    vort2flux!( state.ψ, state.q, state.Γ,
                sfd.ib_prob.model, sfd.ib_prob.model.grid.mg );

    # Store filtered RHS for use in next time step
    for n=length(β):-1:2
        sfd.state.nonlin[n] .= sfd.state.nonlin[n-1]
    end
end
