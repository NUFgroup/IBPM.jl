include("../../src/ibpm.jl")

# Define grid
xlims = (-1.5, 6.5)
ylims = (-2.0, 2.0)
boundary = (xlims..., ylims...) #left, right, bottom, and top of domain
mg = 5   # num domains
Δx = 0.02
grid =  ibpm.make_grid(Δx, boundary, mg=mg)

# Other parameters
Re = 100.0

x0, y0 = 0.0, 0.0  # Plate center
#nb = 0;  # Number of body points (0 for ds=2h)
nb = 51;
L = 1.0   # Plate length

# Initialize motion
Uinf = 1.0;  # Free-stream velocity
α = 0.0;

 #freestream conditions
freestream = (Ux=t->Uinf, Uy=t->0.0, inclination=t->α)

function run_sim!(t, state, prob; output=1, callback=(state, prob)->nothing)
	for i=1:length(t)
		ibpm.advance!(state, prob, t[i])
        if mod(i,output) == 0
			callback(state, prob);  # Primitive callback, can be used for plotting or other output
            @show (t[i], state.CD, state.CL, state.cfl)
        end
	end
end
