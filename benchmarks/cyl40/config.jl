include("../../src/ibpm.jl")

# Define grid
xlims = (-1.0, 3.0)
ylims = (-2.0, 2.0)
boundary = (xlims..., ylims...) #left, right, bottom, and top of domain
mg = 5   # num domains
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
