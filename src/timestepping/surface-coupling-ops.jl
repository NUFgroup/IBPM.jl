"""
    update_surface_quantities!(F̃b, qs, q0, prob)

Solve the Poisson equation (25) in Colonius & Taira (2008).

Dispatch based on the type of motion in the problem - allows precomputing
    regularization and interpolation where possible.
"""
function update_surface_quantities!(qs, state, prob)
    update_surface_quantities!( MotionType(prob.model.bodies),
        BodyType(prob.model.bodies), qs, state, prob)

    return nothing
end

"""
Version for rigid stationary bodies (or rigid bodies with perscribed motion
    simulated in body-fixed frame with appropriate 'fictitious' forces added)

Solve modified Poisson problem for uB = 0 and bc2 = 0
```
 Bf̃ = Eq = ECψ
```
"""
function update_surface_quantities!(::Union{Type{Static}, Type{MovingGrid}},
                        ::Type{RigidBody{T}} where T <: Motion,
                        qs, state, prob)
    Q = prob.model.work.q2  # Net flux

    broadcast!(+, Q, qs, state.q0[:,1])        # qs + q0
    mul!(state.F̃b, prob.model.mats.E, Q)  # E*(qs .+ state.q0)... using fb here as working array
    state.F̃b .= prob.Binv*state.F̃b              # Allocates a small amount of memory

    return nothing
end

"""
Version for rigid bodies with perscribed motion simulated in lab-fixed frame

Solve the Poisson problem for bc2 = 0 (???) with nonzero boundary velocity ub
```
Bf̃ = Eq - ub
   = ECψ - ub
```
"""
function update_surface_quantities!(::Type{V} where V <: Motion,
                          ::Type{RigidBody{T}} where T <: Motion,
                          qs, state, prob)
    # Working memory for in-place operations (small allocation)
    F̃work = similar(state.F̃b)
    Q = prob.model.work.q2   # Net flux

    broadcast!(+, Q, qs, state.q0[:,1])                # qs + q0
    mul!(F̃work, prob.model.mats.E, Q)       # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(state.F̃b, prob.Binv, F̃work)

    return nothing
end

"""
Version for deforming bodies with no additional pprescribed motion (or
    with perscribed motion simulated in body-fixed frame with appropriate
    'fictitious' forces added)

Solve the coupled non-linear system for surface stresses and body displacements
    Iterates via Newton-Raphson to guess consistent stresses and
    displacements/velocities
"""
function update_surface_quantities!(::Union{Type{Static}, Type{MovingGrid}},
                          ::Type{DeformingBody{T}} where T <: Motion,
                          qs, state, prob)

    Q = prob.model.work.q2  # Net flux (aliased into pre-allocated array,
                            # doesn't change over an FSI loop so load here)
                            # needed for getting surface stresses
    broadcast!(+, Q, qs, state.q0[:,1])        # qs + q0

    # FSI loop: N-R to iterate to convergence
    err_FSI = 1.0 #initialize error as large
    tol_FSI = 1.e-5 #tolerance for convergence

    #initialize displacement guesses from values at last time step
    χ_k = copy(state.χ)
    ζ_k = copy(state.ζ)
    ζ̇_k = copy(state.ζ̇)
    f_k = copy(state.fb)

    while err_FSI > tol_FSI

        #Get structure matrices & update E & B matrix
        prob.model.mats.Ms, prob.model.mats.Ks, prob.model.mats.Qs =
                get_structure_mats(typeof(prob.model.bodies[1]),
                                    prob.model.bodies[1])

        #structural matrix needed for time stepping
        #lin solve involves this and B from the flow solver
        K̂s = prob.model.mats.Ks +
            4.0/(prob.scheme.dt^2.0) * prob.model.mats.Ms
        K̂s = inv(K̂s)

        QĨWmat = 2.0 * mats.Qs * mats.Ĩ * mats.W * prob.model.grid.h/prob.scheme.dt^2.0
            #the h and dt scaling comes from the time stepping

        prob.model.mats.E, mats.W = setup_reg( grid, prob.model.bodies )
        prob.Binv = get_Binv(prob.model, prob.Ainv[1], K̂s, QĨWmat)

        #Develop RHS for linearized system for stresses
        mul!(state.F̃b, prob.model.mats.E, Q)  # E*(qs .+ state.q0)... using fb here as working array
        r_c = 2.0/prob.scheme.dt * (χ_k - state.χ) - state.ζ


        r_ζ = prob.model.mats.Ms * ( state.ζ̇ + 4.0/prob.scheme.dt * state.ζ +
            4.0/(prob.scheme.dt^2.0)*(state.χ - χ_k) ) - prob.model.mats.Ks*χ_k

        r_ζ = K̂s * r_ζ

        F_bg = -1.0*prob.model.grid.h*( 2.0/prob.scheme.dt*r_ζ + r_c)
        #This is the RHS in terms of the structural unknowns. Have to
        #re-structure in terms of the fluid surface stress
        F_sm = prob.model.mats.Ĩ * F_bg

        #create the full rhs
        rhsf = F_sm + prob.model.mats.E*Q

        #Solve for surface stresses
        F̃_kp1 = prob.Binv*rhsf
        f_kp1 = prob.model.mats.W*F̃_kp1 / prob.scheme.dt

        #use stressses to get surface disp
        Δχ = r_ζ + K̂s * QĨWmat*F̃_kp1 / prob.scheme.dt

        #update all structural quantities
        χ_k = χ_k + Δχ
        ζ_k = -state.ζ + 2.0/prob.scheme.dt * (χ_k - state.χ)
        ζ̇_k = 4.0/prob.scheme.dt^2.0 * (χ_k - state.χ) - 4.0/prob.scheme.dt *
            state.ζ - state.ζ̇

        #compute error and determine if converged
        err_FSI = norm( f_kp1-f_k )/ norm( f_k )

        #update ds, xb, ...
        f_k = copy(f_kp1)

        for j = 1:length(prob.model.bodies)
            nb, nf = get_body_info([prob.model.bodies[j]])
            prob.model.bodies[j].xb = prob.model.bodies[j].xref +
                reshape(prob.model.mats.Ĩ*χ_k, nb, 2)
            prob.model.bodies[j].ub = reshape(prob.model.mats.Ĩ*ζ_k, nb, 2)
            for jj = 2 : size(xb,1)
                prob.model.bodies[j].ds[jj] = sqrt(
                    (prob.model.bodies[j].xb[jj,1]-prob.model.bodies[j].xb[jj-1,1])^2.0
                    + (prob.model.bodies[j].xb[jj,2]-prob.model.bodies[j].xb[jj-1,2])^2.0 )
            end
            prob.model.bodies[j].ds[1] = prob.model.bodies[j].ds[2]


        end

    end

    for j=1:length(prob.model.bodies)
        state.xb[j] = prob.model.bodies[j].xb
        state.ub[j] = prob.model.bodies[j].ub
    end
end



"""
    update_prescribed_body_motion!(state, prob, t)

Update the immersed bodies and coupling matrices (if applicable).
"""
function update_prescribed_body_motion!(state::IBState, prob::IBProblem, t::Float64)
    update_prescribed_body_motion!(MotionType(prob.model.bodies), state, prob, t)
    return nothing
end

" No motion for static or moving grid cases "
function update_prescribed_body_motion!(::Union{Type{Static}, Type{MovingGrid}},
                        state::IBState,
                        prob::IBProblem,
                        t::Float64)
    return nothing
end

" For moving grid, update the bodies and the coupling operators"
function update_prescribed_body_motion!(::Type{MotionFunction}, state::IBState, prob::IBProblem, t::Float64)
    model = prob.model
    bodies, grid = prob.model.bodies, prob.model.grid
    for j=1:length(bodies)
        move_body!(bodies[j], t)
        state.xb[j] = bodies[j].xb
        state.ub[j] = bodies[j].ub
    end

    # For arbitrary motion in a fixed frame, have to update operators
    model.mats.E, mats.W = setup_reg( grid, bodies )
    prob.Binv = get_Binv(model, prob.Ainv[1])

    return nothing
end
