include("config.jl")  # Set up grid and other parameters

data = load("benchmarks/cyl40_stability/dns_output.jld2")
state = data["state"]

for lev=1:grid.mg
    qx, qy = grid.split_flux(state.q[:, lev])
    qx .= 0.5*(qx .+ qx[:, end:-1:1])
    qy .= 0.5*(qy .- qy[:, end:-1:1])
end
