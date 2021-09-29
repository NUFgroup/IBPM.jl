include("config.jl")  # Set up grid and other parameters

# Load pre-computed DNS solution
using FileIO
data = load("benchmarks/cyl40_stability/dns_output.jld2")
base_state = data["state"]

# Construct the IBProblem corresponding to the DNS
Δt = 1e-2
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

T = 0.1
t = 0:Δt:T
evals, evecs, info = eigsolve(
x₀, 2, :LM, verbosity=3, krylovdim=128, tol=1e-6, orth=ModifiedGramSchmidt()
) do x
    y = deepcopy(x)
    run_sim!(t, y, prob, output=50)
    return y
end

save("benchmarks/cyl40_stability/krylov_output.jld2", "evals", evals, "evecs", evecs, "info", info)
