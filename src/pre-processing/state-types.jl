
"""
State variables (stores everything needed for time stepping)
"""
abstract type State end

mutable struct IBState <: State
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
    IBState() = new() # Default constructor (all undefined references)
end

"Construct empty state from AbstractIBProblem"
function IBState(prob)
    grid = prob.model.grid
    nb, nf = get_body_info(prob.model.bodies)

    state = IBState()
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


"Randomly initialize the vorticity of the IBState"
function IBState(prob, noise_level::Number)
    state = IBState(prob)
    state.Γ = noise_level*randn(size(state.Γ))
    return state
end

function Base.similar(v::IBState)
    w = IBState()
    w.q = similar(v.q)
    w.q0 = similar(v.q0)
    w.Γ = similar(v.Γ)
    w.ψ  = similar(v.ψ)
    w.nonlin = [similar(v.nonlin[i]) for i=1:length(v.nonlin)]
    w.fb = [similar(v.fb[i]) for i=1:length(v.fb)]
    w.F̃b = similar(v.F̃b)
    w.CD = similar(v.CD)
    w.CL = similar(v.CL)
    w.cfl, w.slip = 0.0, 0.0
    return w
end

function Base.:*(α::Number, v::IBState)
    w = similar(v)
    w.Γ .*= α
    w.q .*= α
    w.ψ .*= α
    return w
end

"Out of place scalar multiplication; multiply vector v with scalar α and store the result in w"
function LinearAlgebra.mul!(w::IBState, v::IBState, α::Number)

end
