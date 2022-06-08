```@meta
CurrentModule = LevelSetQuadrature
```

# [LevelSetQuadrature.jl](@id home-section)

*High-order quadratures for implicitly defined surfaces*

!!! note "Algorithm and references"
    This package implements, with minor modifications, the algorithms described
    in:
    - Saye 2015
    - Saye 2022
    If you find this work useful, please consider citing the original
    references. 
    TODO: properly cite.

## Overview

This package allows you to generate high-order *quadrature nodes* and *weights*
for integrating regular functions over $d$-dimensional volumes/surfaces of
the form:

```math
\Omega = \{\boldsymbol{x} \in U \subset \mathbb{R}^d \quad s.t. \quad  \phi(\boldsymbol{x}) < 0\}, \quad \text{and} \quad \Gamma = \{\boldsymbol{x} \in U \subset \mathbb{R}^d \quad s.t. \quad  \phi(\boldsymbol{x}) = 0\},
```

where ``U`` is the bounding box described as a $d$-dimensional rectangle, and
``\phi : \mathbb{R}^d \to \mathbb{R}`` is a *simple Julia function* (more on
this later). We refer to ``\Omega`` and
``\Gamma`` as implicit domains because they are defined *implicitly* through the
function ``\phi``. For example, we may describe a disk through

```@example circle
using LevelSetQuadrature
using Plots

f   = (x) -> x[1]^2 + x[2]^2 - 1
lb, ub = -1.1, 1.1
# plot
xx = yy = lb:0.05:ub
fig = contour(xx,yy,(x,y) -> f((x,y));levels=[0],aspect_ratio=:equal,cbar=nothing)
```

To generate a quadrature for the disk we may then call [`quadgen`](@ref) as follows:

```@example circle
rec = HyperRectangle((lb,lb),(ub,ub))
x,w = quadgen(f,rec,:interior)
println("error in area: ", abs(sum(w) - pi))
```

This produces the following set of nodes

```@example circle
scatter!(fig,PlotPoints(),x,label="interior nodes",ms=3,m=:cross)
```

Similarly, a surface quadrature can be created by passing `:surface` instead of
`:interior` as an argument to [`quadgen`](@ref):

```@example circle
x,w = quadgen(f,rec,:surface)
println("error in perimeter:", abs(sum(w) - 2*pi))
```

```@example circle
scatter!(fig,PlotPoints(),x,label="surface nodes",ms=3)
```

As long as the implicit function defining the domain is simple enough (e.g. a
polynomial), changing the domain should be very straightforward:

```@example cassini
using LevelSetQuadrature
using Plots
# Cassini oval: https://en.wikipedia.org/wiki/Cassini_oval
const a = 1.1
const c = 1.0
f = (x) -> (x[1]^2 + x[2]^2)^2 - 2*c^2*(x[1]^2 - x[2]^2) - (a^4 - c^4)
# a bounding box
l,u = -1.5,1.5
rec = HyperRectangle((l,l),(u,u))
# generate the quadrature
x,w = quadgen(f,rec,:interior)

# plot the domain and quadrature
xx = yy = l:0.05:u
contour(xx,yy,(x,y) -> f((x,y));levels=[0],aspect_ratio=:equal,cbar=nothing,lw=3)
scatter!(PlotPoints(),x,label=nothing,ms=3)
```

There is nothing particularly important about dimension 2 (besides the
convenient fact that *2d* is easy to visualize) -- the algorithms implemented are
dimension agnostic. Computing a surface quadrature for the unit sphere for
instance can be achieved through:

```@example sphere
using LevelSetQuadrature
using Plots
f = (x) -> sum(x.*x) - 1
l,u = -1.5,1.5
rec = HyperRectangle((l,l,l),(u,u,u))
x,w = quadgen(f,rec,:surface)
println("error in surface area: ", sum(w) - 4π)
```

```@example sphere
# plot the nodes
scatter(PlotPoints(),x,label=nothing,ms=3,m=:cross)
```

There are several `keyword` arguments that can be passed to `quadgen` which
provide further control on the generated quadrature; you should consult the
[`quadgen`](@ref) documentation for more details.

!!! warning "Limitation to simple functions"
    In order to build a high-order quadrature, the algorithm implemented in this
    package must build rigorous bounds on the implicit function passed and its
    derivatives. For that purpose, we use the bounded linearization idea of
    Saye, implemented in the [`Linearization`](@ref) type, which performs an
    operation akin to automatic differentiation in order to propagate the bounds
    on the linearization of functions. Because we have only implemented a small
    subset of the "linearization rules", it is possible/likely that the call
    `quadgen(f,rec)` will fail if the algorithm does not manage to propagate the
    bounds on `f` and its gradient; that is, if the call `f(::Linearization)`
    fails. Feel free to open an issue (or a `PR`) if `quadgen` does not work
    with your *simple* function.


## Polynomial domains

In the case of polynomial domains ...

## Advanced usage

## Other limitations and known issues

Here are some know limitations of the package:

- Automatic differentiation
- Singular points
- 

## Index

```@index
```