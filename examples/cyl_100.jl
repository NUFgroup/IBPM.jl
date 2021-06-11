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
body = [ibpm.make_cylinder( r, Δx, 0.0, 0.0; motion=ibpm.Static() )] #If motion is not prescribed then Static() is assumed


#save info
svfc = [ (t,state)->state.CL; (t,state)->state.Γ  ] #vector of functions that
                                                    #specify which data to save
svti = [ [Float64[]]; 0.01 ] #corresponding vector of save instances (1 per function)
                        #each can be prescribed as scalar Float or array
svna = [ :Cl, :γ ]
save_info = (save_fcns = svfc, save_times = svti, save_names = svna)

runtime = IBPM_advance( Re, boundary, body, freestream,
    Δx=Δx, Δt=Δt, T=0.05, plot=false, save_info=save_info )

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
