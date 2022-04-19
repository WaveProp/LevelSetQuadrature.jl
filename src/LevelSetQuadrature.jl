module LevelSetQuadrature

using StaticArrays

include("hyperrectangle.jl")
include("linearization.jl")
include("quadrature.jl")

export HyperRectangle, Linearization, quadratureNodesWeights

end # module
