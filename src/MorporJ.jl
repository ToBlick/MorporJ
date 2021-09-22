module MorporJ

using LinearAlgebra, Convex, COSMO

include("distfunctions.jl")
include("io.jl")
include("barycenters.jl")
include("greedy.jl")

greet() = print("Hello World!")

end # module
