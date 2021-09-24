include("config.jl")  # Set up grid and other parameters

# Load pre-computed DNS solution
using FileIO
data = load("dns_output.jld2")
base_state = data["state"]

# Construct the IBProblem corresponding to the DNS
base_prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

# Construct the linearized problem
prob = ibpm.LinearizedIBProblem(base_state, base_prob, Δt)
base_state, base_prob = nothing, nothing  # Free up memory

ϵ = 1e-6   # Noise level for initial vorticity
state = ibpm.IBState(prob, ϵ);
