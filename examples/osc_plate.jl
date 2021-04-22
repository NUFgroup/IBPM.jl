include("../src/ibpm.jl")
using Plots

# Define grid
nx = 400
ny = 200
mg = 3   # num domains

offx = 4.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 8.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 200.0
Δt = 1e-3

# Initialize motion
motion = ibpm.MovingGrid(t -> sin(t), t -> 0.0, t-> 0.0)

# Create plate
L = 1.0; # Plate length
α = 90.0 * π/180.0;      # Angle of attack of plate
x0, y0 = 0.0, 0.5
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0; motion=motion, n=51 )]

prob = ibpm.IBProblem(grid, bodies, Δt, Re);
state = ibpm.IBState(prob);

T=2.0*2π
t = 0:Δt:T
println(length(t))

# Define a function to plot the plate as a grey line
plot_body(body) = display( plot!(body.xb[:, 1], body.xb[:, 2],
                                color=:grey, lw=5, legend=false) )

# Run simulation and save the animation
anim = ibpm.animated_sim(t, state, prob; output=100) do state, prob
        ibpm.plot_state(state, prob.model.grid, clims=(-5, 5))  # Plot vorticity
        plot_body(prob.model.bodies[1])
end
gif(anim, "examples/osc_plate.gif", fps=30)
