include("../../src/ibpm.jl")


# Define grid
nx = 400
ny = 200
mg = 5   # num domains

offx = 1.5; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 8.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 100.0

x0, y0 = 0.0, 0.0  # Plate center
nb = 0;  # Number of body points (0 for ds=2h)
nb = 51;
L = 1.0   # Plate length

# Initialize motion
Uinf = 1.0;  # Free-stream velocity
α = 0.0      # Free-stream AoA
