include("../src/ibpm.jl")

# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)
len = 4.0  # length of domain in x-direction

# Other parameters
Re = 40.0
Δt = 1e-2

Uinf = 1.0;   # Free-stream flow
r = 0.5; # Cylinder radius
body = (type="cylinder", lengthscale=r, motion="static")
T = 10
runtime = ibpm.IBPM_advance( Re, nx,ny, offx, offy, len,
    mg=mg, body=body, Δt=Δt, Uinf=Uinf, T=T, plot=false )
