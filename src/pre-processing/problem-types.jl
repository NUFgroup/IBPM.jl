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
                       bodies::Array{V, 1} where V <: Body,
                       dt::Float64,
                       Re::Float64;
                       freestream::NamedTuple=(Ux=1.0, Uy=0.0, inclination=0.0)
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
        base_flux!(state, prob, 0.0)  # Initialize base flux at time zero
        return state
    end
end


"""
    base_flux!(state::IBState, prob::IBProblem, t::Float64)
Set background flux based on `prob.model.bodies[].motion`
Assumes same free-stream parameters for all motions (<-- CHANGE THIS)
"""
function base_flux!(state::IBState,
                    prob::IBProblem,
                    t::Float64)
    base_flux!(MotionType(prob.model.bodies), state, prob, t)
end

"Initialize irrotational freestream flux when not time-varying"
function base_flux!(::Type{T} where T <: InertialMotion,
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    Ux = freestream.Ux(t)
    Uy = freestream.Uy(t)
    α = freestream.inclination(t)

    nu = grid.ny*(grid.nx+1);  # Number of x-flux points
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:nu, lev ] .= (Ux*cos(α) - Uy*sin(α))* hc      # x-flux
        state.q0[ nu+1:end, lev ] .= (Ux*sin(α) + Uy*cos(α))*hc  # y-flux
    end
end

"Update time-varying background flux for moving grid"
function base_flux!(::Type{MovingGrid},
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    @assert length(prob.model.bodies) == 1 # Assumes only one body
    grid = prob.model.grid
    motion = prob.model.bodies[1].motion
    XX, YY = prob.model.XX, prob.model.YY;
    nu = grid.ny*(grid.nx+1);  # Number of x-flux points
    nq = grid.nq

    ### Rotational part
    Ω = -motion.θ̇(t)
    α = -motion.θ(t)

    ### Potential flow part (note θ = -α for angle of attack)
    Ux0 = motion.U(t)*cos(α)
    Uy0 = motion.U(t)*sin(α)

    state.q0 .*= 0.0
    for lev=1:grid.mg
        hc = grid.h*2^(lev-1);  # Coarse grid spacing

        ### x-fluxes
        @views state.q0[1:nu, lev] .= YY[:, lev]
        @views state.q0[1:nu, lev] .*= -hc*Ω

        ### y-fluxes
        @views state.q0[(nu+1):nq, lev] .= XX[:, lev]
        @views state.q0[(nu+1):nq, lev] .*= hc*Ω

        ### Irrotational part
        @views state.q0[1:nu, lev] .+= hc*Ux0      # x-flux
        @views state.q0[(nu+1):nq, lev] .+= hc*Uy0  # y-velocity
    end
end
