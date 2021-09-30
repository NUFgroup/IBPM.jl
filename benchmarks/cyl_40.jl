include("../src/ibpm.jl")

# Define grid
xlims = (-1.0, 3.0)
ylims = (-2.0, 2.0)
boundary = (xlims..., ylims...) #left, right, bottom, and top of domain
mg = 4   # num domains
Δx = 0.02
grid =  ibpm.make_grid(Δx, boundary, mg=mg)

# Other parameters
Re = 40.0
Δt = 1e-2

Uinf = 1.0;   # Free-stream flow
r = 0.5; # Cylinder radius

cyls = [ibpm.make_cylinder( r, grid.h, 0.0, 0.0 )]

 #freestream conditions
freestream = (Ux=t->1.0, Uy=t->0.0, inclination=t->0.0)

function run_sim!(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end

prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);

state = ibpm.IBState(prob);
T = 100
t = 0:Δt:T

run_sim!(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed run_sim!(t, state, prob; output=20)
println(runtime)
