include("config.jl")  # Set up grid and other parameters

prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

state = ibpm.IBState(prob);
T = 100
t = 0:Δt:T

run_sim(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed run_sim(t, state, prob; output=20)
println(runtime)

using FileIO
FileIO.save("benchmarks/cyl40/dns_output.jld2",  "state", state)
