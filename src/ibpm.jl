module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using InplaceOps  # @! macro

export IBPM_advance

#Caution, include order matters!
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("pre-processing/pre-processing-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("timestepping/timestepping-include.jl")
include("plotting/plotting-include.jl")

"""
Convenience function to solve the full problem and plot final solution
    for flow over a single cylinder

For more control, just use this as a template - see benchmarks and examples
"""
function IBPM_advance(Re, boundary, body, freestream=(Ux=1.0,);
    Δx=missing, mg=5, Δt=missing, T=20.0*dt, plot=false)

    #--extract user params to sim variables
        Δx, Δt, T, freestream = read_user_vars(Δt, Δx, freestream, Re, T)
    #--

    #--Build flow grid
        grid =  make_grid(Δx, boundary, mg=mg)
    #--

    # #Make body
    # r = body.lengthscale
    #
    # if body.motion == "static"
    #     motion=Static()
    # elseif body.motion == "rot_cyl"
    #     motion=RotatingCyl(body.θ̇)
    # end
    # cyls = [make_cylinder( r, grid.h, 0.0, 0.0; motion=motion )]

    prob = IBProblem(grid, body, Δt, Re, freestream=freestream);
    state = IBState(prob);

	t = 0:Δt:T
    run_sim(t[1:2], state, prob) # Pre-compilation for benchmarking
    runtime = @elapsed run_sim(t, state, prob) #advance to final time

    #plotting
    if plot==true
        plot_state(state, prob.model.grid)
        plot_body(prob.model.bodies[1])
    end

    return runtime
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

function run_sim(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
		if ~all(isfinite.(state.CL))
			break
		end
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
			state.cfl = compute_cfl(state, prob)
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end

function animated_sim(update_plot, t, state, prob;
		nplt=100,
		output=1,
		callback=(state, prob)->nothing)
    n_iter = length(t)÷nplt
    anim = @animate for i=1:n_iter
		sim_idx = (i-1)*nplt.+(1:nplt)
		run_sim( @view(t[sim_idx]), state, prob; output=output, callback=callback )
		if ~all(isfinite.(state.CL))
			break
		end
		update_plot(state, prob)
    end
	return anim
end

end
