include("../src/ibpm.jl")
using Plots

# Define grid
nx = 200
ny = 100
mg = 3   # num domains

offx = 1.5; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 8.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 100.0
Δt = 1e-2

# Initialize motion
Uinf(t) = 1.0;
A = 20.0*π/180  # Final pitch amplitude
tc = 5.0        # Shift of tanh pulse
θ(t) = -A*0.5*(1 + tanh(t-tc))
θ̇(t) = -A*0.5*sech(t-tc)^2;
motion = ibpm.BodyFixed(Uinf, θ, θ̇)

# Create plate
x0, y0 = 0.25, 0.0
nb = 48;  # Number of body points
L = 1.0
α = 0.0
bodies = [ibpm.make_plate( L, α, grid.h, x0, y0, motion )]

prob = ibpm.init_prob(grid, bodies, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux

T=15.0
timesteps = round(Int, T/Δt)
println(timesteps)

function run_sim(it_stop, state, prob)
    nplt = 10
    big_iter = it_stop÷nplt
    anim = @animate for i=1:big_iter
        for j=1:nplt
            it = nplt*i + j
            t = prob.scheme.dt*it
            ibpm.advance!(t, state, prob)
        end
        @show (nplt*i, state.CD, state.CL, state.cfl)
        #ibpm.plot_u(state, prob.model.grid; lev=1)
        ibpm.plot_state(state, prob.model.grid; clims=(-5, 5))
        display( plot!(prob.model.bodies[1].xb[:, 1], prob.model.bodies[1].xb[:, 2],
            color=:grey, lw=5, legend=false) )
    end
    return anim
end

anim = run_sim(timesteps, state, prob) #advance to final time
gif(anim, "pitch_hold.gif", fps = 30)
