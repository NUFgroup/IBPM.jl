include("../src/IBPM.jl")


# Define grid
nx = 400  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
len = 16.0  # length of domain in x-direction
offx = 4.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 4.0; # offset in y dirn (same as offx but in y-dirn)

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

function run_sim(it_stop)
        CL = zeros(it_stop, 3)
        for it=1:it_stop
                t = it*dt
                ibpm.advance!(t, state, prob)
                CL[it, :] = state.CL
                println([it, state.CD, state.CL, state.cfl])
        end
        CL = zeros(it_stop, 3)
end


function run_to_conv(ϵ)
        δCL = 1.0
        CL_prev = zeros(3)
        it = 0
        while δCL > ϵ
                it += 1
                t = it*dt
                CL_prev = copy(state.CL)
                ibpm.advance!(t, state, prob)
                δCL = maximum( abs.(CL_prev .- state.CL) )
                println([it, state.CL[2], δCL, state.cfl])
        end
        return state.CD, state.CL
end

ϵ = 1e-8
CD, CL = run_to_conv(ϵ)
println(CD)
println(CL)
