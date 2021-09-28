include("config.jl")  # Set up grid and other parameters

# Load pre-computed DNS solution
using FileIO
data = load("benchmarks/cyl40/dns_output.jld2")
base_state = data["state"]

# Construct the IBProblem corresponding to the DNS
Δt = 1e-6
base_prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

# Construct the linearized problem
prob = ibpm.LinearizedIBProblem(base_state, base_prob, Δt)
#base_state, base_prob = nothing, nothing  # Free up memory

# using LinearMaps
# L = LinearMap(length(state), ismutating=true) do x
#     ibpm.advance!(y, x, prob, 0.0)
# end
#
# function linearize(base_state, base_prob, Δt; ϵ=1e-6)
#     prob = ibpm.LinearizedIBProblem(base_state, base_prob, Δt)
#     x₀ = ibpm.IBState(prob, ϵ)
#     L = LinearMap(length(state), ismutating=true) do (y, x)
#         ibpm.advance!(y, x, prob, 0.0)
#     end
#     return L, x₀, prob
# end

using KrylovKit
ϵ = 1e-6   # Noise level for initial vorticity
x₀ = ibpm.IBState(prob, ϵ);
sol = eigsolve(x₀, verbosity=3) do x
    y = ibpm.advance(x, prob, 0.0)
    return y
end
