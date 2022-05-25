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
include("bernsteinpolynomials.jl")
include("quadrature.jl")

export
    HyperRectangle,
    Linearization,
    _quadgen,
    quadgen,
    BernsteinPolynomial,
    ∇,
    rebase,
    power2Berstein

end # module
