module LevelSetQuadrature

using LinearAlgebra
using StaticArrays
using Roots
using FastGaussQuadrature

import WavePropBase:
    HyperRectangle,
    section,
    center,
    domain,
    low_corner,
    high_corner,
    half_width,
    UniformCartesianMesh,
    NodeIterator,
    svector

include("linearization.jl")
include("Bernsteinpolynomials.jl")
include("quadrature.jl")

export HyperRectangle, Linearization, quadratureNodesWeights

end # module
