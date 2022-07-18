module LevelSetQuadrature

using LinearAlgebra
using StaticArrays
using Roots
using FastGaussQuadrature
using Requires

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


function __init__()
    @require DynamicPolynomials="7c1d4256-1411-5781-91ec-d7bc3513ac07" begin
        include("dynamicpolynomials_wrapper.jl")
    end
end

end # module
