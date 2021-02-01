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
r = 0.5; # Cyinder radius

# Create an array of one cylinder
motion = ibpm.Static();
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0, motion )];

# Create full IBPM problem
Re = 40.0
dt = 5e-3

χ = 0.47
Δ = 0.33

prob = ibpm.init_prob(grid, cyls, Re, dt);
state = ibpm.init_state(prob);

Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack
ibpm.base_flux!(state, grid, Uinf, α)  # Initialize irrotational base flux

sfd_prob = ibpm.init_sfd(prob, state, Δ, χ)

function run_to_conv(ϵ)
        δCD = 1.0
        CD_prev = 0.
        it = 0
        conv = false
        while (δCD > ϵ)
                it += 1
                t = it*dt
                CD_prev = state.CD[1]
                ibpm.advance!(t, state, sfd_prob)
                δCD = abs.(CD_prev .- state.CD[1])
                println([it, state.CD, state.CL, δCD])
        end
        return state.CD, state.CL
end

ϵ = 1e-8
CD, CL = run_to_conv(ϵ)
