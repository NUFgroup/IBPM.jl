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
Uinf(t) = [sin(t), 0.0];
Ω(t) = 0.0;
motion = ibpm.BodyFixed(Uinf, Ω)

# Create plate
L = 1.0; # Plate length
α = 90.0 * π/180.0;      # Angle of attack
x0, y0 = 0.0, 0.5
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0, motion; n=51 )]

prob = ibpm.init_prob(grid, bodies, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux

T=2.0*2π
timesteps = round(Int, T/Δt)
println(timesteps)


function run_sim(it_stop, state, prob)
    nplt = 50
    big_iter = it_stop÷nplt
    anim = @animate for i=1:big_iter
        for j=1:nplt
            it = nplt*i + j
            t = prob.scheme.dt*it
            ibpm.advance!(t, state, prob)
        end
        @show (nplt*i, state.CD, state.CL, state.cfl)
        ibpm.plot_state(state, prob.model.grid; clims=(-5, 5))
        display( plot!(prob.model.bodies[1].xb[:, 1], prob.model.bodies[1].xb[:, 2],
            color=:grey, lw=10, legend=false) )
    end
    return anim
end

anim = run_sim(timesteps, state, prob) #advance to final time
gif(anim, "osc_plate.gif", fps = 30)
