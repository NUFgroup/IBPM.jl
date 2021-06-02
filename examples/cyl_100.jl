include("../src/ibpm.jl")
using .ibpm

# Problem parameters
boundary = (-1.0453, 3.15, -2.0148, 2.312) #left, right, bottom, and top of domain
Re = 100.0 #Reynolds #
freestream = (Ux=t->t^0.0,) #freestream conditions
                              #can be provided as constants or
                              #functions of time
                              #(default: (Ux=1.0, Uy=0.0, inclination=0.0))

#simulation parameters (these are all optional)
Δx = 0.02 #Default gives a grid Re of 2
Δt = missing #default (==missing) aims for a CFL of 0.1 with a
          #fairly conservative safety factor on max vel
mg=5      #Number of sub-domains. Default is 5



# Create an array of one cylinder
# Create cylinder
r = 0.5; # Cylinder radius
# Create an array of one cylinder
motion = "static"
body = (type="cylinder", lengthscale=r, motion=motion)
runtime = IBPM_advance( Re, boundary, freestream, Δx=Δx,
    body=body, Δt=Δt, T=10.0, plot=true )

# T=300.0
# t = 0:Δt:T
#
# ibpm.run_sim(t[1:2], state, prob) # Pre-compile
#
# # Run simulation and save the animation
# anim = ibpm.animated_sim(t, state, prob; output=20) do state, prob
#         ibpm.plot_state(state, prob.model.grid, clims=(-3, 3))  # Plot vorticity
#         ibpm.plot_cyl(prob.model.bodies[1]);
# end
#
# gif(anim, "examples/cyl_100.gif", fps=10)
