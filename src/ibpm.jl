module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers

export IBPM_advance

#Caution, include order matters!
include("pre-processing/pre-processing-include.jl")
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("plotting/plotting-include.jl")
include("timestepping/timestepping-include.jl")

function IBPM_advance(Re, nx, ny, offx, offy, len; mg=1,body, Δt,
    Uinf=1.0, α=0.0, T=20.0*dt, plot=false)
    # MultiGrid
    grid = make_grid(nx, ny, offx, offy, len, mg=mg)

    #Make body
    r = body.lengthscale

    if body.motion == "static"
        motion=Static()
    elseif body.motion == "rot_cyl"
        motion=RotatingCyl(body.θ̇)
    end
    cyls = [make_cylinder( r, grid.h, 0.0, 0.0, motion )]

    prob = init_prob(grid, cyls, Re, Δt);
    state = init_state(prob);

    println(prob.model.bodies)

    Uinf = []
    base_flux!(state, grid, motion)  # Initialize irrotational base flux

    timesteps = round(Int, T/Δt)

    run_sim(1, state, prob) #pre-compute stationary IB matrix before advancing
    runtime = @elapsed run_sim(timesteps, state, prob) #advance to final time

    #plotting
    if plot==true
        plot_state(state, prob.model.grid)
        plot_body(prob.model.bodies[1])
    end

    return runtime
end


function run_sim(it_stop, state, prob)
    for it=1:it_stop
        t = prob.scheme.dt*it
        advance!(t, state, prob)
        if mod(it,20) == 0
            @show (it, state.CD, state.CL, state.cfl)
        end
    end
end

end
