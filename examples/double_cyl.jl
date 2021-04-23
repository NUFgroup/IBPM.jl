include("../src/ibpm.jl")
using .ibpm

"""
Any[5.0, [2.43850303680305, 2.4758858950482403], [1.8523198640158964, -1.9568651471240626], 0.14081634953573285]
"""

# Define grid
nx = 300  # num of x points on finest domain
ny = 400  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 3   # num domains

offx = 0.8; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.1; # offset in y dirn (same as offx but in y-dirn)
len = 3.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 200.0
Δt = 1e-3

# Initialize motion
Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack

# Create cylinder
r = 0.5; # Cylinder radius
g = 1.0; # Gap between cylinders
# Critical point: (g0 , Re0 ) = (0.725, 56.46)
nb = 0
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, -(g/2+r); n=nb ),
        ibpm.make_cylinder( r, grid.h, 0.0,  (g/2+r); n=nb )]

prob = ibpm.IBProblem(grid, cyls, Δt, Re, Uinf=Uinf);
state = ibpm.IBState(prob);

T=10.0
t = 0:Δt:T

ibpm.run_sim(t[1:2], state, prob) # Pre-compile

# Run simulation and save the animation
anim = ibpm.animated_sim(t, state, prob; output=10, nplt=500) do state, prob
        ibpm.plot_state(state, prob.model.grid, clims=(-5, 5))  # Plot vorticity
        ibpm.plot_cyl(prob.model.bodies[1]);
        ibpm.plot_cyl(prob.model.bodies[2]);
end

gif(anim, "examples/double_cyl.gif", fps=10)
