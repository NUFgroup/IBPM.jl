include("../src/ibpm.jl")
using .ibpm

#Below, the Reynolds #, boundary of the finest domain, and body must be specified.
#all other variables have default values and needn't be provided the user

#--necessary variables
    boundary = (-1.0453, 3.15, -2.0148, 2.312) #left, right, bottom, and top of domain
    Re = 100.0 #Reynolds #

    # specify body as named tuple with keys type, lengthscale, center, motion.
    #type and lengthscale must be specified. The others have defaults associated
    #with a stationary body centered at (x,y)=(0.0,0.0)
    type = :cylinder #:cylinder, :plate are supported
    lengthscale = 0.5 #key lengthscale. e.g., for cylinder is radius. Supports Float64
    motion = :static #type of body motion. supports :static or a function of time
                     #default is static
    center = [0.0; 0.0] #body CoM is centered here. default: [0.0; 0.0]
    body = (type=:cylinder, lengthscale=lengthscale, motion=:static)
#--

#--optional variables
    freestream = (Ux=t->t^0.0,) #freestream conditions
                                #can be provided as constants or
                                #functions of time
                                #(default: (Ux=1.0, Uy=0.0, inclination=0.0))

    T = 0.05 #final time to run to (default = 20.0*dt)

    #simulation parameters (these are all optional)
    Δx = 0.02 #default (==missing) gives a grid Re of 2
    Δt = missing #default (==missing) aims for a CFL of 0.1 with a
              #fairly conservative safety factor on max vel
    mg=5      #Number of sub-domains. Default is 5
#--

#--save_info (optional)
    #gives the code information for a data structure to return
    #user provides as a Named Tuple with three keys: save_fcns, save_times,
    #   save_types
    #default (if save_info is unspecified): the code returns necessary
    #   information for a restart.

    #Even if save_info is provided, only the save_fcns key is necessary.

    #vector of functions that specify which data to save.
    #default: save full state at the final time
    svfc = [ (t,state)->state.CL; (t,state)->state.Γ  ]

    #corresponding vector of save instances (1 per function)
    #each can be prescribed as scalar Float or array.
    #default: save at the final time instance T for each save variable
    #if any of the entries is either 0.0 or missing, that variable will be saved
    #every timestep
    svti = [ 0.0, 0.01 ]

    #types of each variable
    #default: Any for each save variable
    svty = [ Vector{Float64}; Array{Float64,2} ]

    #Store in Named Tuple
    save_info = ( save_fcns = svfc, save_types = svty )

#--

#--run the simulation based on user data
    prob, data, runtime = IBPM_advance( Re, boundary, body, freestream,
        Δx=Δx, Δt=Δt, T=T, plot=false, save_info=save_info )
#--

# T=300.0
# t = 0:Δt:T
#
# ibpm.run_sim(t[1:2], state, prob) # Pre-compile
#
# Run simulation and save the animation
# anim = ibpm.animated_sim(t, state, prob; output=20) do state, prob
#         ibpm.plot_state(state, prob.model.grid, clims=(-3, 3))  # Plot vorticity
#         ibpm.plot_cyl(prob.model.bodies[1]);
# end

# gif(anim, "examples/cyl_100.gif", fps=10)
