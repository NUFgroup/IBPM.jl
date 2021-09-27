include("config.jl")  # Set up grid and other parameters

prob = ibpm.IBProblem(grid, cyls, Δt, Re, freestream=freestream);
state = ibpm.IBState(prob);

Δ = 0.5   # SFD damping
χ = 3.0   # SFD filter width
sfd_prob = ibpm.init_sfd(prob, state, Δ, χ)

T = 100
t = 0:Δt:T

sfd_tol = 1e-8
function sfd_callback(x, sfd)
    x̄ = sfd.state
    ϵ = sqrt( mapreduce(x->x^2, +, x.Γ - x̄.Γ) )
    @show ϵ
    return ϵ < sfd_tol
end

run_sim(t[1:2], state, prob) # Pre-compile

# Advance to final time
runtime = @elapsed run_sim(t, state, sfd_prob; output=20, callback=sfd_callback)
println(runtime)

using FileIO
FileIO.save("sfd_output.jld2",  "state", state)
