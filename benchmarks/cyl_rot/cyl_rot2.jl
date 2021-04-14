"""
Rotating cylinder benchmark 2:

Body-fixed frame

Expected results:
    Converge to CD=2.1562, CL=-0.17185 by T=100
"""

include("../../src/ibpm.jl")
using .ibpm
using MAT

# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 2.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# Other parameters
Re = 20.0
Δt = 1e-2

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Set up motion
Uinf(t) = 1.0;
#θ̇(t) = 0.1;
#θ(t) = θ̇(t)*t;

# For debugging
t_init = 0;  # Let flow equilibrate first
θ̇(t) = t>t_init ? 0.1 : 0.0
θ(t) = θ̇(t)*(t-t_init);
motion = ibpm.MovingGrid(Uinf, θ, θ̇)

# Create cylinder
r = 0.5; # Cylinder radius
cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0, motion; n=78 )]

prob = ibpm.init_prob(grid, cyls, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux

T=100.0
t = 0:Δt:T

function plot_lift(t, F)
	CD = @. F[:, 1]*cos(θ(t)) - F[:, 2]*sin(θ(t))
	CL = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))

	display(plot(t, CL, xlabel="t", ylabel="CL"))#, ylims=(-0.2, 0.0)))
end

function run_sim(t, F, state, prob)
    for i=1:length(t)
        ibpm.advance!(t[i], state, prob)
        F[i, 1] = state.CD[1]
        F[i, 2] = state.CL[1]
        if mod(i,20) == 0
            @show (i, state.CD, state.CL, state.cfl)
			plot_lift(t[1:i], F[1:i, :])
			#ibpm.plot_u(state, grid; lev=3)
			#ibpm.plot_state(state, grid, clims=(-5, 5))
			#ibpm.plot_body(prob.model.bodies[1])
			#CL = F[i, 1]*sin(θ(t[i])) + F[i, 2]*cos(θ(t[i]))
			#println(CL)
        end
    end
end

F = zeros(length(t), 2)
runtime = @elapsed run_sim(t, F, state, prob) #advance to final time

matwrite("benchmarks/cyl_rot/force_jl2.mat", Dict(
	"t" => Array(t),
	"F" => F,
))
