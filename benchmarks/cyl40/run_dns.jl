include("config.jl")  # Set up grid and other parameters

prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

state = ibpm.IBState(prob);
T = 100
t = 0:Δt:T

function run_sim(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end

run_sim(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed run_sim(t, state, prob; output=20)

using FileIO
FileIO.save("dns_output.jld2",  "state", state)
