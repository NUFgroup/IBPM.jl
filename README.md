# IBPM.jl

Julia CFD package based on the Immersed Boundary Projection Method (IBPM) from Taira & Colonius (2007), with multigrid method from Colonius & Taira (2008).

### Overview

The package is generally structured after the hierarchy of the C++ package, with appropriate modifications for Julia (i.e. composite types and multiple dispatch instead of object-oriented).  

* __Bodies__: the bodies are a collection of boundary points at which the no-slip conditions are enforced, along with a Motion (currently only Static "motion" is supported).  The `ib_domain.jl` file also has utilities to make a couple of basic bodies (a cylinder and 4-digit NACA airfoils, for instance).
* __Model__: the Navier-Stokes model is a structure to hold information about the domain - the grid, an array of bodies, and any matrices and operators that can be precomputed.
* __Problem__: Ideally this could interface with the DifferentialEquations.jl package eventually, so a "IBPMProblem" is defined in the same spirit as the "ODEProblem".  This combines a Model with precomputed operators relevant for advancing the state and pre-allocated working memory
* __State__: This is just a collection of pre-allocated arrays to store the vorticity, streamfunction, velocity flux, and body forces (including lift and drag).  It also has storage for the "memory" of the multi-step scheme for the nonlinear terms.  

Once a Problem and State are defined, the basic time-stepping is to just call `advance!(state, prob)`

### Current Status (12/21/2020)
Single- and multi-grid versions have both been tested and optimized.  Currently only fixed, rigid bodies are supported and there is no support for adjoints, automatic differentiation, etc.

I tried parallelizing for loops in `nonlin` and `multigrid_utils` with the `@threads` macro, but it actually ran slower (on my laptop with 4 threads).  I think this is because of shared memory issues.  BLAS should be multi-threaded by default (start julia with `-t {nproc}`)
