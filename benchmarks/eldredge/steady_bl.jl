### Steady flat plate boundary layer
include("setup_prob.jl")  # Set up grid and other variables
using BSON: @save, @load

Δt = 1e-2
T = 100.0
t = 0:Δt:T

# Create plate
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0; n=nb)]
prob = ibpm.init_prob(grid, bodies, Δt, Re; Uinf=Uinf);

#@load "benchmarks/eldredge/steady.bson" state
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, t[1])  # Initialize irrotational base flux

function run_sim(t, state, prob)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
		@show (t[i], state.CD, state.CL, state.cfl)
	end
end

run_sim(t, state, prob) #advance to final time

@save "benchmarks/eldredge/results/steady.bson" state
