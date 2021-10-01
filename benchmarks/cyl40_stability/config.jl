include("../../src/ibpm.jl")

# Define grid
xlims = (-1.0, 7.0)
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

function symmetrize!(state, grid)
	sym(q)  = 0.5*(q .+ q[:, end:-1:1])
	asym(q) = 0.5*(q .- q[:, end:-1:1])
	for lev=1:grid.mg
	    qx, qy = grid.split_flux(@view(state.q[:, lev]))
	    Γ = reshape(@view(state.Γ[:, lev]), grid.nx-1, grid.ny-1)

	    qx .=  sym(qx)
	    qy .= asym(qy)
	    Γ  .= asym(Γ)
	end
end
