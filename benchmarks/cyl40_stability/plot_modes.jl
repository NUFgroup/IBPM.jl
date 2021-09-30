include("config.jl")  # Set up grid and other parameters

x = LinRange(xlims..., grid.nx-1)
y = LinRange(ylims..., grid.ny-1)

using FileIO
# data = load("benchmarks/cyl40_stability/dns_output.jld2")
# base_state = data["state"]
# ω = reshape(base_state.Γ[:, 1], grid.nx-1, grid.ny-1) / grid.h^2
# contourf(x, y, ω', c=cgrad(:seaborn_icefire_gradient),
#         lw=0, levels=10, aspect_ratio=:equal)

data = load("benchmarks/cyl40_stability/krylov_output.jld2")
evals = data["evals"]
evecs = data["evecs"]

ω = reshape(evecs[1].Γ[:, 1], grid.nx-1, grid.ny-1) / grid.h^2
contourf(x, y, real.(ω)', c=cgrad(:seaborn_icefire_gradient),
        lw=0, levels=10, aspect_ratio=:equal)
