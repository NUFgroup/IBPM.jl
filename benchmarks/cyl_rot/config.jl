include("../../src/ibpm.jl")

# Define grid
nx = 200
ny = 200
mg = 5   # num domains
offx = 2.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# MultiGrid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 20.0

x0, y0 = 0.0, 0.0  # Cylinder center
nb = 78;   # Number of body points (0 for ds=2h)
r = 0.5;   # Cylinder radius

# Initialize motion
Uinf = 1.0;  # Free-stream velocity
α = 0.0;

\Omega = 0.1;  # Rotation of the cylinder
