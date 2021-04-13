"""
Rotating cylinder benchmark 3:

Lab frame with moving body: recomputing
    interpolation/regularization operators

Expected results:
"""

include("../src/ibpm.jl")
using .ibpm

# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 2.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# Other parameters
Re = 20.0
Δt = 1e-2
