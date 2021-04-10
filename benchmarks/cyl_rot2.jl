"""
Rotating cylinder benchmark 2:

Body-fixed frame

Expected results:
    Converge to CD=2.1562, CL=-0.17185 by T=100
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

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Set up motion
Uinf(t) = [1.0, 0.0];
Ω(t) = 0.1;
motion = ibpm.BodyFixed(Uinf, Ω)

# Create cylinder
r = 0.5; # Cylinder radius
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0, motion )]

prob = ibpm.init_prob(grid, cyls, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux

T=100.0
timesteps = round(Int, T/Δt)

function run_sim(it_stop, state, prob)
    for it=1:it_stop
        t = prob.scheme.dt*it
        ibpm.advance!(t, state, prob)
        if mod(it,20) == 0
            @show (it, state.CD, state.CL, state.cfl)
        end
    end
end

run_sim(1, state, prob)
runtime = @elapsed run_sim(timesteps, state, prob) #advance to final time
