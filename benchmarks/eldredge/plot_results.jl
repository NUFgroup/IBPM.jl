using Plots
using LaTeXStrings
using MAT
using DelimitedFiles

a = 11
t0 = 0
t1, t2, t3, t4 = 1+t0, 3+t0, 4+t0, 6+t0
α_max = 45*π/180
G(t) = @. log( cosh(a*(t-t1))*cosh(a*(t-t4))/(  cosh(a*(t-t2))*cosh(a*(t-t3)) ) )
Ġ(t) = @. a*(tanh(a*(t-t1))-tanh(a*(t-t2))-tanh(a*(t-t3))+tanh(a*(t-t4)));

file = matopen("benchmarks/eldredge/results/force_jl.mat", "r")
t_jl = read(file, "t")
CL_jl = read(file, "CL")
CD_jl = read(file, "CD")
close(file)

data = readdlm("benchmarks/eldredge/results/force_cpp1.dat")
t_cpp1 = data[:, 2]
#CD_cpp = data[:, 3]
CL_cpp1 = data[:, 4]

# C++ Body-fixed
data = readdlm("benchmarks/eldredge/results/force_cpp2.dat")
t_cpp2 = data[:, 2]
#CD_cpp = data[:, 3]
F = data[:, 3:4]

G_max = maximum(G(t_cpp2))
θ(t) = -α_max*G(t)/G_max
θ̇(t) = -α_max*Ġ(t)/G_max

CL_cpp2 = @. F[:, 1]*sin(θ(t_cpp2)) + F[:, 2]*cos(θ(t_cpp2))
CD_cpp2 = @. F[:, 1]*cos(θ(t_cpp2)) - F[:, 2]*sin(θ(t_cpp2))

data = readdlm("benchmarks/eldredge/results/force_fort.dat")
t_fort = data[:, 1]*1e-3
CD_fort = data[:, 4]
CL_fort = data[:, 5]

plot(t_jl, CL_jl, lw=2, xlims=(0, 7), label="Julia",
    xlabel=L"$t$", ylabel=L"$CL$")
plot!(t_cpp1, CL_cpp1, lw=2, label="C++ (1)")
plot!(t_cpp2, CL_cpp2, lw=2, label="C++ (2)")
plot!(t_fort, CL_fort, lw=2, label="Fortran")
