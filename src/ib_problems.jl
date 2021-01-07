
mutable struct WorkingMemory
    q1::AbstractArray
    q2::AbstractArray
    q3::AbstractArray
    q4::AbstractArray
    Γ1::AbstractArray
    Γ2::AbstractArray
    Γ3::AbstractArray
    bc::AbstractArray
end



struct IBProblem
    model::IBModel
    scheme::ExplicitScheme
    work::WorkingMemory
    A
    Ainv
    B
    Binv
end


function init_memory(grid::UniformGrid)
    return WorkingMemory(
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ)
    )
end

function init_memory(grid::MultiGrid)
    return WorkingMemory(
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ)
    )
end

function init_prob(grid::T where T <: Grid,
                   bodies::Array{V, 1} where V <: Body,
                   Re::Float64,
                   dt::Float64)
    mats = get_mats(grid, bodies, Re)
    model = IBModel(grid, bodies, Re, mats)
    scheme = AB2(dt)   # Explicit time-stepping for nonlinear terms
    A, Ainv, B, Binv = get_AB(model, dt)
    work = init_memory(grid)
    return IBProblem(model, scheme, work, A, Ainv, B, Binv)
end



function MotionType(prob::IBProblem)
    bodies = prob.model.bodies
    motions = [typeof(bodies[i].motion) for i=1:length(bodies)]
    if all(motions .== Static)
        return Static
    else
        return Motion
    end
end
