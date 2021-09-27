include("../../src/ibpm.jl")


# Define grid
nx = 400  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
len = 12.0  # length of domain in x-direction
offx = 3.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 3.0; # offset in y dirn (same as offx but in y-dirn)

# MultiGrid
grid = ibpm.make_grid(nx, ny, offx, offy, len; mg=mg)

# Create cylinder
r = 0.5; # Cylinder radius


# Create an array of one cylinder
motion = ibpm.Static();
x0 = 1.5*√3*r
y0 = 1.5*r
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0, motion),
        ibpm.make_cylinder( r, grid.h,  x0,  y0, motion),
        ibpm.make_cylinder( r, grid.h,  x0, -y0, motion)];

# Create full IBPM problem
Re = 15.0
dt = 5e-3

prob = ibpm.init_prob(grid, cyls, Re, dt);
state = ibpm.init_state(prob);

Uinf = 1.0;   # Free-stream flow
α = 0.0 * π/180.0;      # Angle of attack
ibpm.base_flux!(state, grid, Uinf, α)  # Initialize irrotational base flux

# SFD parameters
χ = 0.47
Δ = 0.33
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
                println([it, state.CD[1], state.CL[1], δCD])
        end
        return state.CD, state.CL
end

ϵ = 1e-8
CD, CL = run_to_conv(ϵ)
println(CD)
println(CL)

ibpm.save_state("output/steady.mat", state);
ibpm.save_model("output/model.mat", prob.model);
