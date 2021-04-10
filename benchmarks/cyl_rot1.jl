"""
Rotating cylinder benchmark 1:

Lab frame, adding "artificial" rotation to boundary points
    i.e. specifying velocity, but not recomputing operators

Expected results:
    Converge to CD=2.1618, CL=-0.211 by T=80
"""

include("../src/ibpm.jl")
using .ibpm

# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# Other parameters
Re = 20.0
Δt = 1e-2
Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack

# Create an array of one cylinder
# Create cylinder
r = 0.5; # Cylinder radius

# Create an array of one cylinder
T=100.0

# MultiGrid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

θ̇ = 0.1;
motion = ibpm.RotatingCyl(Uinf, α, θ̇)
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0, motion )]

prob = ibpm.init_prob(grid, cyls, Re, Δt);
state = ibpm.init_state(prob);

ibpm.base_flux!(state, grid, motion)  # Initialize irrotational base flux
timesteps = round(Int, T/Δt)

ibpm.run_sim(1, state, prob) #pre-compute stationary IB matrix before advancing
runtime = @elapsed ibpm.run_sim(timesteps, state, prob) #advance to final time
