
using LinearAlgebra: norm  # FOR DEBUGGING

#Driver function that steps forward in time
function advance!(t::Float64,
                  state::IBState{UniformGrid},
                  prob::IBProblem)
    grid = prob.model.grid

    # Move bodies and update coupling matrices (E)
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == BodyFixed
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state

    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)
    circ2_st_vflx!( state.ψ, state.q, state.Γ, prob.model );

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    state.cfl = maximum( abs.( (1/grid.h^2) * state.q * dt ) ) ;
end

function advance!(t::Float64,
                  state::IBState{MultiGrid},
                  prob::IBProblem)
    """
    After updating get_nonlin:
        84.443 ms (1961 allocations: 76.93 MiB)
    After optimization:
        40.325 ms (719 allocations: 15.66 MiB)
    """
    grid = prob.model.grid
    update_bodies!(prob, t)

    if MotionType(prob.model.bodies) == BodyFixed
        base_flux!(state, prob, t)
    end

    # Alias working memory for notational clarity
    #   This leaves work.Γ2, Γ3 and work.q1, q2 available
    qs = prob.work.q1  # Trial flux
    Γs = prob.work.Γ1  # Trial circulation

    #Computes trial circulation Γs and associated strmfcn and vel flux that
    #don't satisfy no-slip (from explicitly treated terms)
    get_trial_state!(qs, Γs, state, prob)

    # Update surface quantities to be able to trim off part of circ
    # that doesn't satisfy no slip
    @views boundary_forces!(state.F̃b, qs[:, 1], state.q0[:, 1], prob)
    update_stress!(state, prob) #Compute integral quantities and store in state

    # --update circulation , vel-flux, and strmfcn on fine grid
    #   to satisfy no-slip updates state.Γ, state.ψ, state.q
    project_circ!(Γs, state, prob)    # OK TO HERE

    # --Update circulation on all grids based on fine-grid correction
    for lev = 2:grid.mg
        @views coarsify!( state.Γ[:, lev-1], state.Γ[:, lev], grid);
    end

    circ2_st_vflx!( state.ψ, state.q, state.Γ, prob.model, grid.mg);

    #--A few simulation quantities of interest
    # get CFL (u * dt / dx) :
    dt = prob.scheme.dt
    state.cfl = maximum( abs.( (1/(grid.h^2)) * @view(state.q[:, 1]) * dt ) ) ;
end

function base_flux!(state::IBState,
                    prob::IBProblem,
                    t::Float64)
    base_flux!(MotionType(prob.model.bodies), state, prob, t)
end

function base_flux!(::Type{T} where T <: InertialMotion,
                    state::IBState{UniformGrid},
                    prob::IBProblem,
                    t::Float64)
    """
    Initialize irrotational freestream flux

    Assumes same free-stream parameters for all motions
    """
    grid, motion = prob.model.grid, prob.model.bodies[1].motion
    Uinf, α = motion.Uinf, motion.α
    m = grid.nx;
    n = grid.ny;
    state.q0[ 1:(m-1)*n ] .= Uinf * grid.h * cos(α);  # x-flux
    state.q0[ (m-1)*n+1:end ] .= Uinf * grid.h * sin(α);  # y-flux
end

function base_flux!(::Type{T} where T <: InertialMotion,
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    """
    Initialize irrotational freestream flux

    Assumes same free-stream parameters for all motions
    """
    grid, motion = prob.model.grid, prob.model.bodies[1].motion
    Uinf, α = motion.Uinf, motion.α
    m = grid.nx;
    n = grid.ny;
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:(m-1)*n, lev ] .= Uinf * hc * cos(α);      # x-flux
        state.q0[ (m-1)*n+1:end, lev ] .= Uinf * hc * sin(α);  # y-flux
    end
end

function base_flux!(::Type{BodyFixed},
                    state::IBState{UniformGrid},
                    prob::IBProblem,
                    t::Float64)
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

function base_flux!(::Type{BodyFixed},
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    motion = prob.model.bodies[1].motion
    nx, ny, h = grid.nx, grid.ny, grid.h;
    nu = ny*(nx-1);  # Number of x-flux points

    ### Rotational part
    Ω = -motion.θ̇(t)
    α = -motion.θ(t)

    state.q0 .*= 0.0
    for lev=1:grid.mg
        hc = h*2^(lev-1);  # Coarse grid spacing

        ### x-flux
        y = @. ((1:ny)-0.5-ny/2)*hc + ny/2*h - grid.offy
        YY = ones(nx-1)*y'
        state.q0[1:nu, lev] .= -hc*Ω*YY[:]

        ### y-flux
        x = @. ((1:nx)-0.5-nx/2)*hc + nx/2*h - grid.offx
        XX = x*ones(ny-1)'
        state.q0[nu+1:end, lev] .= hc*Ω*XX[:]
    end

    ### Potential flow part (note θ = -α for angle of attack)
    Ux0 = motion.U(t)*cos(α)
    Uy0 = motion.U(t)*sin(α)
    for lev=1:grid.mg
        # Coarse grid spacing
        hc = h*2^(lev-1 );

        state.q0[1:nu, lev] .+= hc*Ux0          # x-flux
        state.q0[nu+1:end, lev] .+= hc*Uy0  # y-velocity
    end
end

"""
#TODO: DO THIS WITHOUT circ2_st_vflx (DIRECT CROSS PRODUCT)
function base_flux!(::Type{BodyFixed},
                    state::IBState{MultiGrid},
                    prob::IBProblem,
                    t::Float64)
    grid = prob.model.grid
    motion = prob.model.bodies[1].motion
    nu = grid.ny*(grid.nx-1);  # Number of x-flux points

    # Alias working memory
    ψ0 = prob.work.Γ1;
    Γ0 = prob.work.Γ2;

    ### Rotational part
    # Get flux from curl of solid-body rotation
    #  Note factor of 2 from angular velocity -> vorticity
    Ω = -motion.θ̇(t)
    α = -motion.θ(t)

    for lev=1:grid.mg
        hc = grid.h * 2^( lev - 1 );  # Coarse grid spacing
        Γ0[:, lev] .= 2*hc^2*Ω
    end

    circ2_st_vflx!(ψ0, state.q0, Γ0, prob.model, grid.mg)

    ### Potential flow part (note θ = -α for angle of attack)
    Ux0 = motion.U(t)*cos(α)
    Uy0 = motion.U(t)*sin(α)
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        ### Potential flow part
        state.q0[1:nu, lev] .+= hc*Ux0          # x-flux
        state.q0[nu+1:grid.nq, lev] .+= hc*Uy0  # y-velocity
    end
end"""

function get_trial_state!(qs::AbstractArray,
                        Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    Combine explicit Laplacian and nonlinear terms into a rhs
       then use Ainv to return trial circulation Γs

    High-level version of AB2:
    rhs = A*Γ .-
    #    3*dt/2 * nonlin .+
    #    dt/2 * nonlin_prev

    Then do Ainv of that to back out trial circ
    """
    dt = prob.scheme.dt
    work = prob.work
    rhs = work.Γ2  # RHS of discretized equation

    #compute the nonlinear term for the current time step
    get_nonlin!( state.nonlin[1], state, prob );

    # Explicit part of Laplacian
    mul!(rhs, prob.A, state.Γ)

    # Explicit nonlinear terms from multistep scheme
    for n=1:length(prob.scheme.β)
        work.Γ3 .= state.nonlin[n]
        rmul!(work.Γ3, prob.scheme.β[n]*dt)
        rhs .-= work.Γ3
    end

    # Store current nonlinear term
    state.nonlin[2] .= state.nonlin[1];

    # Trial circulation  Γs = Ainv * rhs
    mul!(Γs, prob.Ainv, rhs);

    # Trial velocity  (note ψ is used here as a dummy variable)
    circ2_st_vflx!( state.ψ, qs, Γs, prob.model );
end

function get_trial_state!(qs::AbstractArray,
                        Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    Combine explicit Laplacian and nonlinear terms into a rhs
       then use Ainv to return trial circulation Γs

    High-level version of AB2:
    rhs = A*Γ .-
    #    3*dt/2 * nonlin .+
    #    dt/2 * nonlin_prev .+
    #    dt/2 * rhsbc

    Then do Ainv of that to back out trial circ
    """
    dt = prob.scheme.dt
    work = prob.work
    grid = prob.model.grid
    rhsbc = work.bc
    rhs = work.Γ2  # RHS of discretized equation

    for lev = grid.mg:-1:1

        # compute the nonlinear term for the current time step
        get_nonlin!( @view(state.nonlin[1][:, lev]), state, prob, lev );

        # contribution of Laplacian term...
        rhsbc .*= 0.0
        if lev < grid.mg
            # from explicit treatment of circulation
            get_Lap_BCs!( rhsbc, @view(state.Γ[:, lev+1]), lev, prob.model )

            # from current (trial) vorticity at previous grid level
            get_Lap_BCs!( rhsbc, @view(Γs[:,lev+1]), lev, prob.model)
        end

        # Combine explicit Laplacian and nonlinear terms into a rhs
        @views mul!(rhs[:, lev], prob.A[lev], state.Γ[:, lev])

        for n=1:length(prob.scheme.β)
            work.Γ3[:, lev] .= state.nonlin[n][:, lev]
            rmul!(work.Γ3, prob.scheme.β[n]*dt)
            rhs .-= work.Γ3
        end

        # Include boundary conditions
        #   High-level: rhs += 0.5*rhsbc
        work.Γ3[:, lev] .= rhsbc
        rmul!(work.Γ3, 0.5*dt)
        rhs .+= work.Γ3

        # Trial circulation  Γs = Ainv * rhs
        # TODO: use @view to do in-place multiplication
        # Doesn't work here because of view and FFT plan for even indices... WHY??
        Γs[:, lev] .= prob.Ainv[lev] * rhs[:, lev]
    end

    # Store nonlinear solution for use in next time step
    state.nonlin[2] .= state.nonlin[1]

    # Trial velocity on 1st grid level (don't need all grids)
    #   Streamfunction state.ψ is a dummy variable here to compute qs from Γs
    grid = prob.model.grid
    circ2_st_vflx!( state.ψ, qs, Γs, prob.model, 2 );
end

"""
Solve the modified Poisson equation (26)
Dispatch based on the type of motion in the problem
Allows precomputing regularization and interpolation for Static motions
"""
function boundary_forces!(F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    boundary_forces!( MotionType(prob.model.bodies), F̃b, qs, q0, prob)
end

"""
Solve the Poisson equation (25) in Colonius & Taira (2008)
    for uB = 0 and bc2 = 0
    Bf̃ = Eq
       = ECψ
"""
function boundary_forces!(::Union{Type{Static}, Type{BodyFixed}},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    E = prob.model.mats.E
    h = prob.model.grid.h
    # Working memory for in-place operations
    qwork = @view(prob.work.q2[:, 1])
    broadcast!(+, qwork, qs, q0)                     # qs + q0
    mul!(F̃b, E, qwork)                               # E*(qs .+ state.q0)... using fb here as working array
    F̃b .= (1/h)*prob.Binv*F̃b                         # Allocates a small amount of memory
end

function boundary_forces!(::Type{RotatingCyl},
                          F̃b::AbstractVector,
                          qs::AbstractVector,
                          q0::AbstractVector,
                          prob::AbstractIBProblem)
    """
    Solve the Poisson equation (25) in Colonius & Taira (2008)
        for bc2 = 0 with specialized situation of rotating cylinder
    In this case the points don't need to move, but do have nonzero velocity
        Bf̃ = Eq - ub
           = ECψ - ub
    """
    E = prob.model.mats.E
    h = prob.model.grid.h

    # Working memory for in-place operations
    F̃work = similar(F̃b)
    qwork = @view(prob.work.q2[:, 1])

    broadcast!(+, qwork, qs, q0)           # qs + q0
    mul!(F̃work, E, qwork)                  # E*(qs .+ state.q0)
    F̃work .-= get_ub(prob.model.bodies)*prob.model.grid.h   # Enforce no-slip conditions
    mul!(F̃b, prob.Binv, F̃work/h);
end



function project_circ!(Γs::AbstractArray,
                     state::IBState,
                     prob::IBProblem)
    """
    Update circulation to satisfy no-slip condition

    Dispatch based on the type of motion in the problem
    Allows precomputing regularization and interpolation for Static motions
    """
    project_circ!(MotionType(prob.model.bodies), Γs, state, prob)
end

function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{UniformGrid},
                       prob::IBProblem)
    """
    High-level version:
        state.Γ[:, 1] .= Γs .- Array(prob.Ainv[1] * (mats.RET*fb_til_dt) )
    """
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    RET = prob.model.mats.RET   # Precomputed R * E'
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs   # Now Γs is free for working memory
    @views mul!( Γs[:, 1], RET, fb_til_dt)
    @views mul!( Γwork, prob.Ainv, Γs[:, 1])
    state.Γ[:, 1] .-= Γwork
end


function project_circ!(::Type{V} where V<:Motion,
                       Γs::AbstractArray,
                       state::IBState{MultiGrid},
                       prob::IBProblem)
    """
    High-level version:s
        state.Γ[:, 1] .= Γs .- Array(prob.Ainv[1] * (mats.RET*fb_til_dt) )
    """
    Γwork = @view(prob.work.Γ3[:, 1]) # Working memory
    RET = prob.model.mats.RET   # Precomputed R * E'
    fb_til_dt = state.F̃b

    # Low-level version:
    state.Γ .= Γs
    @views mul!( Γs[:, 1], RET, fb_til_dt)
    @views mul!( Γwork, prob.Ainv[1], Γs[:, 1])  # This is the only difference with the UniformGrid version
    state.Γ[:, 1] .-= Γwork
end

#Utilities for storing stress values
function update_stress!(state::IBState,
                        prob::IBProblem)
    """
    Store surface stresses and integrated forces

    Mutates "state"
    """
    nb, nf = get_body_info(prob.model.bodies)
    h = prob.model.grid.h
    dt = prob.scheme.dt
    fb_til_dt = state.F̃b

    # Store surface stress and integrated forces
    nbod_tally = 0; #  Used to keep a tally of which body we're on
    for j = 1 : length(prob.model.bodies)
        ds = prob.model.bodies[j].ds
        # surface stresses
        state.fb[j] .= fb_til_dt[nbod_tally .+ (1:nf[j])] *(h / dt) ./ [ds; ds] ;

        # integrated forces
        state.CD[j] = 2 * sum( ds .* state.fb[j][1 : nb[j]] ) ;
        state.CL[j] = 2 * sum( ds .* state.fb[j][1 .+ nb[j] : nf[j] ] ) ;

        # update body index
        nbod_tally += nf[j];
    end
end


function AB2(dt::Float64)
    """
    Initialize second-order Adams-Bashforth scheme
    """
    return AdamsBashforth(dt, [1.5, -0.5])
end

#TODO: break out by multiple dispatch... but don't duplicate code
function update_bodies!(prob::IBProblem, t::Float64)
    model = prob.model
    bodies, grid = prob.model.bodies, prob.model.grid
    motion = MotionType(bodies)

    if motion != Static
        for j=1:length(bodies)
            move_body!(bodies[j], t)
        end
    end

    # For arbitrary motion in an inertial frame, have to update operators
    if motion == MotionFunction
        model.mats.E = coupling_mat( grid, bodies )
        model.mats.RET = (model.mats.E*model.mats.C)'
        prob.Binv = get_B(model, prob.Ainv)
    end
end
