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
                       freestream::NamedTuple=(Ux=t->0.0,Uy=t->0.0,inclination=t->0.0*t^0.0)
                       )
        prob = new()
        prob.model = IBModel(grid, bodies, Re; freestream=freestream)
        prob.scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        prob.A, prob.Ainv, prob.Binv = get_AB(prob.model, dt)
        return prob
    end
end

"""
Modified IBProblem to include base state.  Only modification to the code
is the direct product called by the `nonlinear!` function
"""
mutable struct LinearizedIBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    base_state::IBState
    QB::Array{Float64, 2}
    ΓB::Array{Float64, 2}
    A
    Ainv
    Binv
    """
        LinearizedIBProblem(base_state, base_prob)

    Create a linearized problem from the base_state and associated IBProblem

    NOTE: The `freestream` value in the nonlinear IBProblem will be incorrect
        for the linearized case, but this field is not used in
        base_flux!(..., prob::LinearizedIBProblem, ...)
    """
    function LinearizedIBProblem(
                    base_state::IBState,
                    base_prob::IBProblem,
                    dt::Float64)
        prob = new()

        prob.model = deepcopy(base_prob.model)
        prob.scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
        prob.A, prob.Ainv, prob.Binv = get_AB(prob.model, dt)
        prob.base_state = deepcopy(base_state)

        # Averaged base flux used in mean flow advection
        #  Note this calls the IBProblem version of avg_flux, not Linearized
        #  so that the "background" or free-stream flux is accounted for
        prob.QB = zeros(prob.model.grid.nq, prob.model.grid.mg)
        for lev=1:prob.model.grid.mg
            prob.QB[:, lev] = copy(avg_flux(base_state, base_prob; lev=lev))
        end
        prob.ΓB = copy(base_state.Γ)

        return prob
    end
end
