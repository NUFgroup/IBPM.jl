"""
IBProblem has the info of IBModel as well as the problem structure (e.g., the
explicit time stepping scheme and information about the implicit treatment via
the A and B matrices and their inverses).

Maybe some opportunity for restructuring...

Looking towards possible compatibility with DifferentialEquations.jl, this
        would be similar to the ODEProblem
"""

abstract type AbstractIBProblem end

"""
Initialize the problem structure (matrices used, bodies and simulation
parameters, time steppping scheme, ...)

Note: the scheme actually speaks to the terms that are explicitly treated. This
is a projection method the directly enforces the no-slip condition, so some terms
are implicitly treated. This information is not contained in scheme, but in the
A, Ainv, B, and Binv matrices
"""
mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    A
    Ainv
    Binv
    function IBProblem(grid::T where T <: Grid,
                       bodies::Array{<:Body, 1},
                       dt::Float64,
                       Re::Float64;
                       freestream::NamedTuple
                       )
        prob = new()
        prob.model = IBModel(grid, bodies, Re; freestream=freestream)
        prob.scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        prob.A, prob.Ainv, prob.Binv = get_AB(prob.model, dt)
        return prob
    end
end


"""
State variables (stores everything needed for time stepping)
"""
abstract type State end

mutable struct IBState{T<:Grid} <: State
    q::Array{Float64, 2}
    q0::Array{Float64, 2}
    Γ::Array{Float64, 2}     # Circulation
    ψ::Array{Float64, 2}     # Streamfunction
    nonlin::Array{Array{Float64, 2}, 1}  # Memory of nonlinear terms
    fb::Array{Array{Float64, 2}, 1}          # Surface stresses
    F̃b::Array{Float64, 2}                    # Body forces * dt
    CD::Array{Float64, 1}    # Drag coefficient
    CL::Array{Float64, 1}    # Lift coefficient
    cfl::Float64
    slip::Float64
    xb::Array{Array{Float64, 2}, 1}
    function IBState(prob::IBProblem)
        grid = prob.model.grid
        nb, nf = get_body_info(prob.model.bodies)

        state = new{typeof(grid)}()
        state.q  = zeros(grid.nq, grid.mg)    # Flux
        state.q0 = zeros(grid.nq, grid.mg)    # Background flux
        state.Γ  = zeros(grid.nΓ, grid.mg)    # Circulation
        state.ψ  = zeros(grid.nΓ, grid.mg)    # Streamfunction
        state.nonlin = [zeros(grid.nΓ, grid.mg) for i=1:length(prob.scheme.β)]
        state.fb = [zeros(nb[i], 2) for i=1:length(nf)]
        state.F̃b = zeros(sum(nf), 1)
        state.CD = zeros(length(prob.model.bodies))
        state.CL = zeros(length(prob.model.bodies))
        state.cfl, state.slip = 0.0, 0.0

        state.xb = [zeros(nb[i], 2) for i=1:length(nf)]
        for j=1:length(prob.model.bodies)
            state.xb[j] = prob.model.bodies[j].xb
        end
        base_flux!(state, prob, 0.0)  # Initialize base flux at time zero
        return state
    end
end
