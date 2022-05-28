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
    svector,
    order

include("linearization.jl")
include("bernsteinpolynomials.jl")
include("quadrature.jl")

export
    HyperRectangle,
    quadgen,
    # integrate,
    BernsteinPolynomial,
    power2bernstein

end # module
