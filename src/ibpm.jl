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
    Δx=missing, mg=5, Δt=missing, T=20.0*dt, plot=false, save_info=missing)

    #--extract user params to sim variables
        Δx, Δt, T, freestream = read_user_vars(Δt, Δx, freestream, Re, T)
    #--

    #--Build flow grid
        grid =  make_grid(Δx, boundary, mg=mg)
    #--

	#--Initialize problem types and models
	    prob = IBProblem(grid, body, Δt, Re, freestream=freestream);
	    state = IBState(prob);
	#--

	#Time over which simulation will be run
	t = 0:Δt:T

	#initialize user desired save variables as a NamedTuple called data
	data = init_save_data( t, save_info, state )

	#run simulation over desired time window
    run_sim(t[1:2], state, prob) # Pre-compilation for benchmarking
    runtime = @elapsed run_sim(t, state, prob, data=data) #advance to final time

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

function run_sim(t, state, prob;
	display_freq=20,
	data::Union{NamedTuple,Vector{Float64}}=Float64[])
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
		if ~all(isfinite.(state.CL))
			break
		end
        if mod(i,display_freq) == 0
			state.cfl = compute_cfl(state, prob)
            @show (t[i], state.CD, state.CL, state.cfl)
        end

		#save stuff if desired
	    save_data!( t[i], t, prob, state, data  )
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
