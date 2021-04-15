using Plots
using MAT
using DelimitedFiles

θ̇(t) = 0.1
θ(t) = @. θ̇(t)*t

file = matopen("benchmarks/cyl_rot/force_jl2.mat", "r")
t_jl = read(file, "t")
F = read(file, "F")
close(file)

CD_jl = @. F[:, 1]*cos(θ(t)) - F[:, 2]*sin(θ(t))
CL_jl = @. F[:, 1]*sin(θ(t)) + F[:, 2]*cos(θ(t))

data = readdlm("benchmarks/cyl_rot/force_cpp1.dat")
t_cpp1 = data[:, 2]
CD_cpp1 = data[:, 3]
CL_cpp1 = data[:, 4]

plot(t_cpp1, CL_cpp1)
plot!(t_jl, CL_jl)
