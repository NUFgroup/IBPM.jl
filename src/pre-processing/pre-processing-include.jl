"""
Information about which files to include within this directory
"""

#include("grid.jl")    # Moved to fluid-domain
#include("motions.jl")  # Moved to structure-domain
#include("bodies.jl")   # Moved to structure-domain
include("models.jl")          # model, matrices, working memory
include("problem-types.jl")   # problem, state, and base flux
include("read-user-vars.jl")
