module LevelSetQuadrature

using StaticArrays
using Plots

include("hyperrectangle.jl")
include("linearization.jl")
include("quadrature.jl")

export HyperRectangle, Linearization, quadratureNodesWeights

end # module
