"""
Information about which files to include within this directory
"""

#include("grid.jl")    # Moved to fluid-domain
include("motions.jl")  # Move to structure-domain?
include("bodies.jl")   # Move to structure-domain?
include("models.jl")          # model, matrices, working memory
include("problem-types.jl")   # problem and state
