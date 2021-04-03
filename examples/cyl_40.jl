include("../src/ibpm.jl")
using .ibpm

# Problem parameters
boundary = (-1.0453, 3.15, -2.0148, 2.312) #left, right, bottom, and top of domain
Re = 40.0 #Reynolds #
freestream = (Ux=t->t^0.0,) #freestream conditions
                              #can be provided as constants or
                              #functions of time
                              #(default: (Ux=1.0, Uy=0.0, inclination=0.0))

#simulation parameters (these are all optional)
Δx = 0.02 #Default gives a grid Re of 2
Δt = missing #default aims for a CFL of 0.1 with a
          #fairly conservative safety factory on max vel
mg=5      #Number of sub-domains. Default is 5


# Create an array of one cylinder
# Create cylinder
r = 0.5; # Cylinder radius
# Create an array of one cylinder
motion = "static"
body = (type="cylinder", lengthscale=r, motion=motion)
runtime = IBPM_advance( Re, boundary, freestream, Δx=Δx,
    body=body, Δt=Δt, T=1.0, plot=true )
