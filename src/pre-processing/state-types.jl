
"""
State variables (stores everything needed for time stepping)
"""
abstract type State end

mutable struct IBState{T<:Number} <: State
    q::Array{T, 2}
    q0::Array{T, 2}
    Γ::Array{T, 2}     # Circulation
    ψ::Array{T, 2}     # Streamfunction
    nonlin::Array{Array{T, 2}, 1}  # Memory of nonlinear terms
    fb::Array{Array{Float64, 2}, 1}          # Surface stresses
    F̃b::Array{Float64, 2}                    # Body forces * dt
    CD::Array{Float64, 1}    # Drag coefficient
    CL::Array{Float64, 1}    # Lift coefficient
    cfl::Float64
    slip::Float64
    xb::Array{Array{Float64, 2}, 1}
    IBState() = new{Float64}() # Default constructor (all undefined references)
    IBState(::Type{T}) where T <: Number = new{T}() # Default constructor (all undefined references)
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
    state.q = noise_level*randn(size(state.q))
    return state
end

"Define the length of a state as the size of the circulation vector Γ"
Base.length(v::IBState) = size(v.Γ, 1)

"Copy all fields of IBState v to w"
function Base.copy!(w::IBState{T}, v::IBState{T}) where T
    # Matrices can be copied normally
    for name in [:q, :q0, :Γ, :ψ, :F̃b, :CD, :CL]
        copyto!( getfield(w, name), getfield(v, name) )
    end
    # Arrays of arrays have to be copied by each element
    for name in [:nonlin, :fb, :xb]
        vfield, wfield = getfield(v, name), getfield(w, name)
        for i=1:length(vfield)
            copyto!(wfield[i], vfield[i])
        end
    end
    w.cfl, w.slip = v.cfl, v.slip
    return w
end

function Base.similar(v::IBState{T}) where T
    return similar(T, v)
end

function Base.similar(::Type{T}, v::IBState{V}) where {T<:Number, V<:Number}
    # Will convert the field type if necessary (saves ops if not)
    sim_field = (T==V) ? field -> similar(field) : field -> convert.(T, similar(field))

    w = IBState(T)  # All fields are :undef, can't be accessed with getfield
    w.q = sim_field(v.q)
    w.q0 = sim_field(v.q0)
    w.Γ = sim_field(v.Γ)
    w.ψ  = sim_field(v.ψ)
    w.nonlin = [sim_field(v.nonlin[i]) for i=1:length(v.nonlin)]
    w.fb = [similar(v.fb[i]) for i=1:length(v.fb)]
    w.F̃b = similar(v.F̃b)
    w.CD = similar(v.CD)
    w.CL = similar(v.CL)
    w.xb = copy(v.xb)
    w.cfl, w.slip = 0.0, 0.0
    return w
end

"Out of place scalar multiplication; multiply vector v with scalar α and store the result in w"
function LinearAlgebra.mul!(w::IBState, v::IBState, α::Number)
    broadcast!(x -> x*α, w.Γ, v.Γ)
    broadcast!(x -> x*α, w.q, v.q)
    broadcast!(x -> x*α, w.ψ, v.ψ)
    for i=1:length(v.nonlin)
        broadcast!(x -> x*α, w.nonlin[i], v.nonlin[i])
    end
    return w
end

"In-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector"
function LinearAlgebra.rmul!(v::IBState, α::Number)
    v.Γ .*= α
    v.q .*= α
    v.ψ .*= α
    for i=1:length(v.nonlin)
        v.nonlin[i] .*= α
    end
    return v
end

"multiply v with a scalar α, which can be of a different scalar type"
function Base.:*(α::T, v::IBState{V}) where {T <: Number, V <: Number}
    w = similar(T, v)
    mul!(w, v, α)
    return w
end
Base.:*(v::IBState{V}, α::T) where {T <: Number, V <: Number} = Base.:*(α, v)

"store in w the result of α*v + β*w"
function LinearAlgebra.axpby!(α::Number, v::IBState, β::Number, w::IBState)
    rmul!(w, β/α)

    w.q .+= v.q
    w.Γ .+= v.Γ
    w.ψ .+= v.ψ
    for i=1:length(v.nonlin)
        w.nonlin[i] .+= v.nonlin[i]
    end

    rmul!(w, α)
    return w
end

"store in w the result of α*v + w"
function LinearAlgebra.axpy!(α::Number, v::IBState, w::IBState)
    LinearAlgebra.axpby!(α, v, one(α), w)
    return w
end

"compute the inner product (in kinetic energy) of two states"
LinearAlgebra.dot(v::IBState, w::IBState) = dot(v.q[:, 1], w.q[:, 1])

"compute the 2-norm (kinetic energy) of a state"
LinearAlgebra.norm(v::IBState) = norm(v.q[:, 1])
