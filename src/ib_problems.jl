abstract type AbstractIBProblem end

mutable struct IBProblem <: AbstractIBProblem
    model::IBModel
    scheme::ExplicitScheme
    work::WorkingMemory
    A
    Ainv
    Binv
end

function init_prob(grid::T where T <: Grid,
                   bodies::Array{V, 1} where V <: Body,
                   Re::Float64,
                   dt::Float64)
    mats = get_mats(grid, bodies, Re)
    model = IBModel(grid, bodies, Re, mats)
    scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
    work = init_memory(grid)
    A, Ainv, Binv = get_AB(model, dt)
    return IBProblem(model, scheme, work, A, Ainv, Binv)
end
