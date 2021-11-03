using JuMP, MathProgBase, CPLEX, LinearAlgebra
include("JuMPeR_master/src/JuMPeR.jl") # change the path accordingly

include("inputStructs.jl")
include("stochasticProgFcts.jl")
include("newsVendorAffinePolicy.jl")
include("newsVendorLiftedPolicy.jl")

include("testCase.jl")