module LevelSetQuadrature

using LinearAlgebra
using StaticArrays
using Roots
using FastGaussQuadrature

import WavePropBase:
    HyperRectangle,
    section,
    center,
    width,
    domain,
    low_corner,
    high_corner,
    half_width,
    UniformCartesianMesh,
    NodeIterator,
    svector,
    order,
    ambient_dimension,
    split,
    PlotPoints

include("linearization.jl")
include("bernsteinpolynomials.jl")
include("implicitdomain.jl")
include("quadrature.jl")

export
    HyperRectangle,
    quadgen,
    # integrate,
    BernsteinPolynomial,
    power2bernstein,
    # re-export
    WavePropBase,
    PlotPoints

end # module
