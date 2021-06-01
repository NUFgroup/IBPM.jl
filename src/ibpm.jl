module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using JLD2

export IBPM_advance

#Caution, include order matters!
include("pre-processing/pre-processing-include.jl")
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("plotting/plotting-include.jl")
include("timestepping/timestepping-include.jl")

function IBPM_advance(Re, boundary, freestream=(Ux=1.0,);
    Δx=missing, mg=5,body, Δt=missing, T=20.0*dt, plot=false)

    #--extract user params to sim variables
        Δx, Δt, T, freestream = read_user_vars(Δt, Δx, freestream, Re, T)
    #--

    #--Build flow grid
        # MultiGrid
        grid =  make_grid(Δx, boundary, mg=mg)
    #--

    #--Build body
        r = body.lengthscale

        if body.motion == "static"
            motion=Static()
        elseif body.motion == "rot_cyl"
            motion=RotatingCyl(body.θ̇)
        end
        cyls = [make_cylinder( r, grid.h, 0.0, 0.0, motion )]
    #--

    #--Initialize problem and state definitions
        prob = init_prob(grid, cyls, Re, Δt, freestream)
        state = init_state(prob)
    #--

    #--Timestepping
        timesteps = round(Int, T/Δt)

        run_sim(1, state, prob) #pre-compute IB matrix before advancing
        runtime = @elapsed run_sim(timesteps, state, prob) #march to final time

        @save "./save_vars.jld2" {compress=true} state prob
    #--

    #--plotting
        if plot==true
            plot_state(state, prob.model.grid)
            for j = 1:length(prob.model.bodies)
                plot_body(prob.model.bodies[1])
            end
        end
    #--

    return runtime
end


function run_sim(it_stop, state, prob)
    for it=1:it_stop
        t = prob.scheme.dt*it
        advance!(t, state, prob)
        if mod(it,20) == 0
            @show (it, state.CD[it,1], state.CL, state.cfl)
        end
    end
end

end
