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
    work::WorkingMemory
    A
    Ainv
    Binv
    function IBProblem(grid::T where T <: Grid,
                       bodies::Array{V, 1} where V <: Body,
                       dt::Float64,
                       Re::Float64;
                       Uinf::Float64=1.0,
                       α::Float64=0.0)
        prob = new()
        prob.model = IBModel(grid, bodies, Re; Uinf=Uinf, α=α)
        prob.scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        prob.A, prob.Ainv, prob.Binv = get_AB(prob.model, dt)
        prob.work = WorkingMemory(grid)
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
    fb::Array{Array{Float64, 1}, 1}          # Surface stresses
    F̃b::Array{Float64, 1}                    # Body forces * dt
    CD::Array{Float64, 1}    # Drag coefficient
    CL::Array{Float64, 1}    # Lift coefficient
    cfl::Float64
    slip::Float64
    function IBState(prob::IBProblem)
        grid = prob.model.grid
        mg = (grid isa UniformGrid) ? 1 : grid.mg  # Number of grid levels
        nb, nf = get_body_info(prob.model.bodies)

        state = new{typeof(grid)}()
        state.q  = zeros(grid.nq, mg)    # Flux
        state.q0 = zeros(grid.nq, mg)    # Background flux
        state.Γ  = zeros(grid.nΓ, mg)    # Circulation
        state.ψ  = zeros(grid.nΓ, mg)    # Streamfunction
        state.nonlin = [zeros(grid.nΓ, mg) for i=1:length(prob.scheme.β)]
        state.fb = [zeros(nf[i]) for i=1:length(nf)]
        state.F̃b = zeros(sum(nf))
        state.CD = zeros(length(prob.model.bodies))
        state.CL = zeros(length(prob.model.bodies))
        state.cfl, state.slip = 0.0, 0.0
        return state
    end
end
