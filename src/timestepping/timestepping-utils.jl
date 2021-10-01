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
    Ux = prob.model.freestream.Ux(t)
    Uy = prob.model.freestream.Uy(t)
    α = prob.model.freestream.inclination(t)

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
    Ux0 = motion.U(t)*cos(α) - motion.V(t)*sin(α)
    Uy0 = motion.U(t)*sin(α) + motion.V(t)*cos(α)

    ## Add in underlying freestream components
    Uxf = prob.model.freestream.Ux(t)
    Uyf = prob.model.freestream.Uy(t)
    αf = prob.model.freestream.inclination(t)

    Ux0 += Uxf*cos(αf)-Uyf*sin(αf)
    Uy0 += Uxf*sin(αf)+Uyf*cos(αf)

    state.q0 .*= 0.0
    for lev=1:grid.mg

        hc = grid.h*2.0^(Float64(lev)-1.0)  # Coarse grid spacing

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


"""
    update_stress!(state, prob)

Store surface stresses and integrated forces.

Mutates "state"
"""
function update_stress!(state::IBState,
                        prob::IBProblem)
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    state.F̃b = prob.model.mats.W * state.F̃b
    state.F̃b = reshape(state.F̃b, sum(nb), 2)

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j = 1 : length(prob.model.bodies)

        ds = prob.model.bodies[j].ds

        # surface stresses
        state.fb[j] .= state.F̃b[nbod_tally .+ (1:nb[j]), :]  / dt

        # integrated forces
        state.CD[j] = 2.0 * sum( ds .* state.fb[j][:, 1] ) ;
        state.CL[j] = 2.0 * sum( ds .* state.fb[j][:, 2] ) ;

        # update body index
        nbod_tally += nb[j];
    end

    state.F̃b = reshape(state.F̃b, 2*sum(nb), 1)
    return nothing
end


"""
compute_cfl(state, prob)

Compute the CFL number (uΔt/Δx) based on the fine-grid flux

Note that this uses working memory that is also used in `nonlinear!`
"""
function compute_cfl(state, prob)
	Δt, Δx = prob.scheme.dt, prob.model.grid.h
	qwork = prob.model.work.q5
	@views @. qwork = abs( state.q[:, 1] )
	return maximum(qwork)*Δt/Δx
end
