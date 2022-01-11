module MorporJ

using LinearAlgebra, Convex, COSMO, Printf, Interpolations

include("distfunctions.jl")
include("io.jl")
include("barycenters.jl")
include("greedy.jl")
include("interpolations.jl")

end # module
