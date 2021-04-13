include("../src/ibpm.jl")
using Plots

# Define grid
nx = 400
ny = 200
mg = 3   # num domains

offx = 1.5; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 8.0  # length of domain in x-direction

# Initialize grid
grid = ibpm.make_grid(nx, ny, offx, offy, len, mg=mg)

# Other parameters
Re = 200.0
Δt = 1e-3

# Initialize motion
Uinf(t) = 1.0;
fb = 0.1
A = 40.0 * π/180.0  # Pitch amplitude, degrees
θ(t) = -A*sin(2π*fb*t)
θ̇(t) = -2π*fb*A*cos(2π*fb*t);
motion = ibpm.BodyFixed(Uinf, θ, θ̇)

# Create plate
x0 = 0.25
nb = 48;  # Number of body points
spec = "0012"
bodies = [ibpm.make_naca(x0, nb, spec, motion)]

prob = ibpm.init_prob(grid, bodies, Re, Δt);
state = ibpm.init_state(prob);
ibpm.base_flux!(state, prob, 0.0)  # Initialize irrotational base flux

T=2.0/fb
timesteps = round(Int, T/Δt)
println(timesteps)

function plot_naca(body)
    nb = size(body.xb, 1)÷2
    xC = body.xb[1:nb, 1]
    yU = body.xb[1:nb, 2]
    yL = body.xb[end:-1:nb+1, 2]
    yC = 0.5*(yU + yL)  # Camber line
    yT = 0.5*(yU - yL)  # Thickness
    plot!(xC, yC, ribbon=yT,
        color=:grey, lw=0, fillalpha=1.0, legend=false)
end

function run_sim(it_stop, state, prob)
    nplt = 100
    big_iter = it_stop÷nplt
    anim = @animate for i=1:big_iter
        for j=1:nplt
            it = nplt*i + j
            t = prob.scheme.dt*it
            ibpm.advance!(t, state, prob)
        end
        @show (nplt*i, state.CD, state.CL, state.cfl)
        ibpm.plot_state(state, prob.model.grid; clims=(-5, 5))
        display( plot_naca(prob.model.bodies[1]) )
    end
    return anim
end

anim = run_sim(timesteps, state, prob) #advance to final time
gif(anim, "pitching_naca.gif", fps = 30)
