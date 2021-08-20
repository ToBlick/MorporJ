module MorporJ

using LinearAlgebra, Convex, COSMO

include("distfunctions.jl")
include("io.jl")
include("barycenters.jl")

greet() = print("Hello World!")

end # module
