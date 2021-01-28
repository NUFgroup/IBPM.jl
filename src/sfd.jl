"""
Selective frequency damping

Implemented as a wrapper around the standard time-stepper
"""

# Only need circulation and a time history of rhs term
#    for the auxiliary filtered state
mutable struct SFDState{T<:Grid} <: State
        Γ::Array{Float64, 2}     # Circulation
        nonlin::Array{Array{Float64, 2}, 1}  # Memory of nonlinear terms
end


mutable struct SFDProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    work::WorkingMemory
    A
    Ainv
    Binv
    Δ::Float    # Filter width
    χ::Float    # Damping
    sfd::State  # Filtered state
end


function init_sfd(grid::T where T <: Grid,
                  bodies::Array{V, 1} where V <: Body,
                  Re::Float64,
                  dt::Float64,
                  x::State,
                  Δ::Float64,
                  χ::Float64)
    prob = init_prob(grid, bodies, Re, dt)
    x̄ = SFDState(copy(x.Γ), copy(x.nonlin))
    return SFDProblem(prob.model,
                      prob.scheme,
                      prob.work,
                      prob.A,
                      prob.Ainv,
                      prob.Binv,
                      x̄, Δ, χ)
end


function get_nonlin!( nonlin::AbstractVector,
                      state::IBState,
                      sfd::SFDProblem,
                      lev::Int )

    # First call the standard routine
    get_nonlin!(nonlin, state, sfd.ib_prob, lev)

    # Will save this in the auxiliary state
    sfd_rhs = @view(sfd.state.nonlin[1][:, lev])

    # Compute filter term χ*(Γ - Γ̄)
    @views broadcast!(-, sfd_rhs, state.Γ[:, lev], sfd.state.Γ[:, lev])
    rmul!(sfd_rhs, sfd.χ)

    # Update nonlinear term
    nonlin .-= sfd_rhs

    # Finally, rescale for updating auxiliary state
    rmul!(sfd_rhs, 1/(sfd.Δ*sfd.χ))
end


function advance_sfd!(state::IBState,
                      sfd::SFDProblem)
      scheme = sfd.ib_prob.scheme
      rhs = sfd.ib_prob.work.Γ3

      # Explicit time stepping for auxiliary state
      for n=1:length(scheme.β)
          rhs .= sfd.state.nonlin[n]
          rmul!(rhs, scheme.β[n]*dt)
          sfd.state.Γ .-= rhs
      end

      # Store nonlinear solution for use in next time step
      for n=1:length(scheme.β)-1
          sfd.state.nonlin[n+1] .= sfd.state.nonlin[n]
      end
  end

function advance!(t::Float64,
                  state::IBState,
                  sfd::SFDProblem)
    # Update primary state
    advance!(t, state, sfd)

    # Update auxiliary state
    advance_sfd!(state, sfd)
end
