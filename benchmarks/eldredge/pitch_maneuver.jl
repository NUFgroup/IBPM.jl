### Benchmark with canonical Eldredge pitch/plunge maneuver
include("config.jl")  # Set up grid and other parameters
using Plots
using BSON: @save, @load
using MAT

Δt = 1e-3
T=10.0
t = 0:Δt:T

### Initialize motion
a = 11
t1, t2, t3, t4 = 1, 3, 4, 6
α_max = 45*π/180
G(t) = @. log( cosh(a*(t-t1))*cosh(a*(t-t4))/(  cosh(a*(t-t2))*cosh(a*(t-t3)) ) )
Ġ(t) = @. a*(tanh(a*(t-t1))-tanh(a*(t-t2))-tanh(a*(t-t3))+tanh(a*(t-t4)));
G_max = maximum(G(t))
θ(t) = -α_max*G(t)/G_max
θ̇(t) = -α_max*Ġ(t)/G_max

motion = ibpm.MovingGrid(t -> Uinf, θ, θ̇)

# Create plate
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0; motion=motion, n=nb)]
prob = ibpm.IBProblem(grid, bodies, Δt, Re, Uinf=Uinf);

# Load steady boundary layer solution from plate_steady.jl
@load "benchmarks/eldredge/results/steady.bson" state

function plot_lift(t, F)
	CD = @. F[:, 1]*cos(θ(t)) - F[:, 2]*sin(θ(t))
	CL = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))

	display(plot(t, CL, xlabel="t", ylabel="CL"))
end

function run_sim(t, F, state, prob)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
		F[i, 1] = state.CD[1]
		F[i, 2] = state.CL[1]

		CD = @. F[i, 1]*cos(θ(t[i])) - F[i, 2]*sin(θ(t[i]))
		CL = @. F[i, 1]*sin(θ(t[i])) + F[i, 2]*cos(θ(t[i]))
		@show (t[i], CD, CL, state.cfl)
	end
end

function animated_sim(t, F, state, prob)
    nplt = 100
    n_iter = length(t)÷nplt
    anim = @animate for i=1:n_iter
		sim_idx = (i-1)*nplt.+(1:nplt)
		run_sim( @view(t[sim_idx]), @view(F[sim_idx, :]), state, prob )
		ibpm.plot_state(state, grid, clims=(-5, 5))
        display( plot!(prob.model.bodies[1].xb[:, 1], prob.model.bodies[1].xb[:, 2],
            color=:grey, lw=5, legend=false) )
    end
	return anim
end

# Hold the body-fixed frame forces
F = zeros(length(t), 2)

# To benchmark
runtime = @elapsed run_sim(t, F, state, prob)

# To visualize
#anim =  animated_sim(t, F, state, prob) #advance to final time
#gif(anim, "benchmarks/eldredge/results/eldredge.gif", fps = 30)

CD = @. F[:, 1]*cos(θ(t)) - F[:, 2]*sin(θ(t))
CL = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))
matwrite("benchmarks/eldredge/results/force_jl2.mat", Dict(
	"t" => Array(t),
	"CL" => CL,
	"CD" => CD,
))
