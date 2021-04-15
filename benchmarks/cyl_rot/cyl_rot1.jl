"""
Rotating cylinder benchmark 1:

Lab frame, adding "artificial" rotation to boundary points
    i.e. specifying velocity, but not recomputing operators

Expected results:
"""
include("config.jl")  # Set up grid and common variables

# Other parameters
T=100.0
Δt = 1e-2

motion = ibpm.RotatingCyl(\Omega)
cyls = [ibpm.make_cylinder( r, grid.h, x0, y0; motion=motion, n=nb )]

prob = ibpm.init_prob(grid, cyls, Δt, Re, Uinf=Uinf);
state = ibpm.init_state(prob);

ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux
timesteps = round(Int, T/Δt)

ibpm.run_sim(1, state, prob) #pre-compute stationary IB matrix before advancing
runtime = @elapsed ibpm.run_sim(timesteps, state, prob) #advance to final time
