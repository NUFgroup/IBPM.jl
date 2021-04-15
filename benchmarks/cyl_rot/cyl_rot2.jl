"""
Rotating cylinder benchmark 2:

Body-fixed frame

Expected results:
"""
using MAT
include("config.jl")  # Set up grid and common variables

motion = ibpm.MovingGrid(t -> Uinf,
						 t -> \Omega*t,
						 t -> \Omega)

# Create cylinder
cyls = [ibpm.make_cylinder( r, grid.h, x0, y0; motion=motion, n=nb )]

prob = ibpm.init_prob(grid, cyls, Δt, Re);
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
