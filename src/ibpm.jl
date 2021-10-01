module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers
using InplaceOps  # @! macro
using Plots

export IBPM_advance

#Caution, include order matters!
include("fluid-domain/fluid-domain-include.jl")
include("structure-domain/structure-domain-include.jl")
include("interface-coupling/interface-coupling-include.jl")
include("pre-processing/pre-processing-include.jl")
include("fluid-operators/fluid-operators-include.jl")
include("structure-operators/structure-operators-include.jl")
include("timestepping/timestepping-include.jl")
include("plotting/plotting-include.jl")

"""
Convenience function to solve the full problem and plot final solution
    for flow over a single cylinder

For more control, just use this as a template - see benchmarks and examples
"""
function IBPM_advance(Re, boundary, body, freestream=(Ux=t->0.0,);
    Δx=missing, mg=5, Δt=missing, T=20.0*dt, save_info=missing)

    #--extract user params to sim variables
        Δx, Δt, T, freestream = read_user_vars(Δt, Δx, freestream, Re, T)
    #--

    #--Build flow grid
        grid =  make_grid(Δx, boundary, mg=mg)
    #--

    #-build body
        body = make_body( body, Δx )
    #--

	#--Initialize problem types and models
	    prob = IBProblem(grid, body, Δt, Re, freestream=freestream)
	    state = IBState(prob)
	#--

	#Time over which simulation will be run
	t = 0:Δt:T

	#initialize user desired save variables as a NamedTuple called data
	data = init_save_data( t, save_info, state )

	#run simulation over desired time window
    runtime = @elapsed run_sim!(t, state, prob, data=data) #advance to final time

    return prob, data, runtime
end


function run_sim!(t, state, prob;
	display_freq=25,
	data::Array{user_var, 1})
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

end #module
