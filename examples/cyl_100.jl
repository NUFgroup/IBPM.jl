include("../src/ibpm.jl")
using .ibpm

# Define grid
# NOTE: This should be approximately doubled for high-quality results
nx = 225  # num of x points on finest domain
ny = 100  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 3   # num domains

offx = 2.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 9.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 100.0
Δt = 1e-2

# Initialize motion
Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack

# Create cylinder
r = 0.5; # Cylinder radius
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0 )]

prob = ibpm.IBProblem(grid, cyls, Δt, Re, Uinf=Uinf);
state = ibpm.IBState(prob);

T=300.0
t = 0:Δt:T

ibpm.run_sim(t[1:2], state, prob) # Pre-compile

# Run simulation and save the animation
anim = ibpm.animated_sim(t, state, prob; output=20) do state, prob
        ibpm.plot_state(state, prob.model.grid, clims=(-3, 3))  # Plot vorticity
        ibpm.plot_cyl(prob.model.bodies[1]);
end

gif(anim, "examples/cyl_100.gif", fps=10)
