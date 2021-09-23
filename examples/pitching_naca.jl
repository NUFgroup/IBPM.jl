include("../src/ibpm.jl")
using Plots

xlims = (-1.5, 6.5)
ylims = (-2.0, 2.0)
boundary = (xlims..., ylims...) #left, right, bottom, and top of domain
mg = 3   # num domains
Δx = 0.02
grid =  ibpm.make_grid(Δx, boundary, mg=mg)

# # Define grid
# nx = 400
# ny = 200
# mg = 3   # num domains

# offx = 1.5; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
# offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
# len = 8.0  # length of domain in x-direction

# # Initialize grid
# grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 200.0
Δt = 1e-3

# Initialize motion
ω = 2π*0.1          # Pitch frequency
A = 40.0 * π/180.0  # Pitch amplitude, degrees
θ(t) = -A*sin(ω*t)
θ̇(t) = -ω*A*cos(ω*t);
U(t) = 1.0
V(t) = 0.0
motion = ibpm.MovingGrid(U, V, θ, θ̇)

# Create plate
x0 = 0.25
nb = 48;  # Number of body points
spec = "0012"
bodies = [ibpm.make_naca(x0, nb, spec, motion=motion)]

prob = ibpm.IBProblem(grid, bodies, Δt, Re);
state = ibpm.IBState(prob);

T=2.0*(2π/ω)
t=0:Δt:T

function plot_naca(body)
    nb = size(body.xb, 1)÷2
    xC = body.xb[1:nb, 1]
    yU = body.xb[1:nb, 2]
    yL = body.xb[end:-1:nb+1, 2]
    yC = 0.5*(yU + yL)  # Camber line
    yT = 0.5*(yU - yL)  # Thickness
    display( plot!(xC, yC, ribbon=yT,
        color=:grey, lw=0, fillalpha=1.0, legend=false) )
end

function run_sim(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end

function animated_sim(update_plot, t, state, prob;
    nplt=100,
    output=1,
    callback=(state, prob)->nothing)
    n_iter = length(t)÷nplt
    anim = @animate for i=1:n_iter
        sim_idx = (i-1)*nplt.+(1:nplt)
        run_sim( @view(t[sim_idx]), state, prob; output=output, callback=callback )
        update_plot(state, prob)
    end
    return anim
end

anim = animated_sim(t, state, prob; output=100) do state, prob
    ibpm.plot_state(prob, state, t, var=:omega, xlims=xlims, ylims=ylims, clims=(-5, 5))  # Plot vorticity
    plot_naca(prob.model.bodies[1])
end
gif(anim, "examples/pitching_naca.gif", fps = 30)
