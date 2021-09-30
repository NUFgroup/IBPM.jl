include("../src/ibpm.jl")

# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 1   # num domains
offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 4.0  # length of domain in x-direction

# Other parameters
Re = 40.0
Δt = 1e-2

Uinf = 1.0;   # Free-stream flow
r = 0.5; # Cylinder radius

grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0 )]

prob = ibpm.IBProblem(grid, cyls, Δt, Re, Uinf=Uinf, α=0.0);
state = ibpm.IBState(prob);
T = 100
t = 0:Δt:T

ibpm.run_sim(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed ibpm.run_sim(t, state, prob; output=20)
