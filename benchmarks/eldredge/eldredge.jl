### Benchmark with canonical Eldredge pitch/plunge maneuver
include("../../src/ibpm.jl")
using Plots

# Define grid
nx = 300
ny = 200
mg = 5   # num domains

offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 6.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 100.0
Δt = 1e-3

T=16.0
t = 0:Δt:T

# Initialize motion
Uinf(t) = 1.0;

a = 11
t0 = 5
t1, t2, t3, t4 = 1+t0, 3+t0, 4+t0, 6+t0
α_max = 45*π/180
G(t) = @. log( cosh(a*(t-t1))*cosh(a*(t-t4))/(  cosh(a*(t-t2))*cosh(a*(t-t3)) ) )
Ġ(t) = @. a*(tanh(a*(t-t1))-tanh(a*(t-t2))-tanh(a*(t-t3))+tanh(a*(t-t4)));
G_max = maximum(G(t))
θ(t) = -α_max*G(t)/G_max
θ̇(t) = -α_max*Ġ(t)/G_max

motion = ibpm.BodyFixed(Uinf, θ, θ̇)

# Create plate
x0, y0 = 0.0, 0.0
nb = 48;  # Number of body points
L = 1.0
α = 0.0
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0, motion; n=nb)]

prob = ibpm.init_prob(grid, bodies, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux


function plot_lift(t, F)
	CD = @. F[:, 1]*cos(θ(t)) - F[:, 2]*sin(θ(t))
	CL = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))

	display(plot(t, CL, xlabel="t", ylabel="CL"))
end

function run_sim(t, F, state, prob)
    nplt = 100
    big_iter = length(t)÷nplt
	iter = 0
    anim = @animate for i=1:big_iter
		for j=1:nplt
			iter += 1
        	ibpm.advance!(t[iter], state, prob)
        	F[iter, 1] = state.CD[1]
        	F[iter, 2] = state.CL[1]
		end
		CL = @. F[iter, 1]*sin(θ(t[iter])) + F[iter, 2]*cos(θ(t[iter]))
        @show (t[iter], state.CD, CL, state.cfl)
		#plot_lift(t[1:i], F[1:i, :])
		ibpm.plot_state(state, grid, clims=(-5, 5))
        display( plot!(prob.model.bodies[1].xb[:, 1], prob.model.bodies[1].xb[:, 2],
            color=:grey, lw=5, legend=false) )
    end
	return anim
end

F = zeros(length(t), 2)
anim =  run_sim(t, F, state, prob) #advance to final time
gif(anim, "eldredge.gif", fps = 30)


# NOTE: Possible difference in force calculation between C++ and Fortran
"""
# Rescaling:
CL = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))
h, ds = grid.h, prob.model.bodes[1].ds[1]
plot(t, CL*(ds/(2*h)), xlims=(5, 12))
"""
