include("../src/IBPM.jl")


# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# MultiGrid
grid = ibpm.make_grid(nx, ny, offx, offy, len; mg=mg)

# Create cylinder
r = 0.5; # Cylinder radius

# Create an array of one cylinder
motion = ibpm.Static();
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, motion )];

# Create full IBPM problem
Re = 40.0
dt = 5e-3

prob = ibpm.init_prob(grid, cyls, Re, dt);
state = ibpm.init_state(prob);

Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack
ibpm.base_flux!(state, grid, Uinf, α)  # Initialize irrotational base flux

function run_sim(it_stop)
        for it=1:it_stop
                ibpm.advance!(state, prob)
                println([it, state.CD, state.CL, state.cfl])
        end
end


run_sim(1)  # First step to compile
println(@elapsed run_sim(100))
