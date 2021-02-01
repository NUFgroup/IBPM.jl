module ibpm

using LinearAlgebra
using SparseArrays
using FFTW
using LinearMaps
using IterativeSolvers

import Base.Threads.@threads

load_plots = true

# Optimized DST for implicit part of Laplacian
include("dst_inversion.jl")

#include("naca.jl")         # For 4-digit NACA airfoils
include("ib_domain.jl")     # Grid, with helpful utilities
include("ib_motions.jl")    # Motion types for bodies
include("ib_bodies.jl")     # Immersed body types, with example constructors

include("ib_matutils.jl")   # Helpful functions for indexing matrices
include("ib_coupling.jl")   # Regularization and interpolation functions
include("ib_mats.jl")       # Precompute various operators

# The "model" is a combination of the grid, bodies, and Reynolds number
#   Initializing the model precomputes matrices from ib_mats.jl
include("models.jl")

# Utilities specific to multigrid formulation (e.g. coarsify)
include("multigrid_utils.jl")

include("schemes.jl")  # Time-stepping schemes for nonlinear term
include("timestepper_utils.jl")  # Utilities for "A" and "B" matrices... things that depend on dt

# The "problem" combines the model with a time-stepping scheme and pre-alllocated memory
#   in the vein of "Problems" from the DifferentialEquations.jl package
include("ib_problems.jl")

# The "state" includes circulation, flux, and streamfunction, along with memory
#    for the history of nonlinear terms (for multistep schemes)
include("ib_state.jl")

include("nonlin.jl")         # Functions to compute the nonlinear terms
include("ibpm_solver.jl")    # Functions to actually advance the state

include("sfd.jl")   # Selective frequency damping

if load_plots
    include("plot_utils.jl")
end


end
