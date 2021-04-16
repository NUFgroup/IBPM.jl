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
                    state::IBState{UniformGrid},
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    Uinf, α = prob.model.Uinf, prob.model.α
    m = grid.nx;
    n = grid.ny;
    state.q0[ 1:(m+1)*n ] .= Uinf * grid.h * cos(α);  # x-flux
    state.q0[ (m+1)*n+1:end ] .= Uinf * grid.h * sin(α);  # y-flux
end

"Initialize irrotational freestream flux when not time-varying"
function base_flux!(::Type{T} where T <: InertialMotion,
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    Uinf, α = prob.model.Uinf, prob.model.α
    m = grid.nx;
    n = grid.ny;
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:(m+1)*n, lev ] .= Uinf * hc * cos(α);      # x-flux
        state.q0[ (m+1)*n+1:end, lev ] .= Uinf * hc * sin(α);  # y-flux
    end
end

"Update time-varying background flux for moving grid"
function base_flux!(::Type{MovingGrid},
                    state::IBState{UniformGrid},
                    prob::IBProblem,
                    t::Float64)
    @assert length(prob.model.bodies) == 1 # Assumes only one body
    grid = prob.model.grid
    motion = prob.model.bodies[1].motion
    nu = grid.ny*(grid.nx-1);  # Number of x-flux points

    ### Rotational part
    Ω = -motion.θ̇(t)
    α = -motion.θ(t)
    nx, ny, h = grid.nx, grid.ny, grid.h;

    ### x-flux
    # TODO: CHECK OFFSETS FOR THESE
    y = h*(1:ny) .- grid.offy
    YY = ones(nx-1)*y'
    state.q0[1:nu, 1] .= -h*Ω*YY[:]

    ### y-flux
    x = h*(1:nx) .- grid.offx
    XX = x*ones(ny-1)'
    state.q0[nu+1:grid.nq, 1] .= h*Ω*XX[:]

    ### Potential flow part
    Ux0 = motion.U(t)*cos(α)
    Uy0 = motion.U(t)*sin(α)
    state.q0[1:nu, 1] .+= h*Ux0          # x-flux
    state.q0[nu+1:grid.nq, 1] .+= h*Uy0  # y-velocity
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
    nu = grid.ny*(grid.nx-1);  # Number of x-flux points

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
        #y = @. ((1:ny)-0.5-ny/2)*hc + ny/2*h - grid.offy
        #YY = ones(nx-1)*y'   # TODO:  Pre-compute or find a better way to do this
        state.q0[1:nu, lev] .= -hc*Ω*YY[:, lev]

        ### y-fluxes
        #x = @. ((1:nx)-0.5-nx/2)*hc + nx/2*h - grid.offx
        #XX = x*ones(ny-1)'   # TODO:  Pre-compute or find a better way to do this
        state.q0[nu+1:end, lev] .= hc*Ω*XX[:, lev]

        ### Irrotational part
        state.q0[1:nu, lev] .+= hc*Ux0          # x-flux
        state.q0[nu+1:end, lev] .+= hc*Uy0  # y-velocity
    end
end
