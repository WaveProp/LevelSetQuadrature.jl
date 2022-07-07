var documenterSearchIndex = {"docs":
[{"location":"references/#references-section","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Modules = [LevelSetQuadrature]","category":"page"},{"location":"references/#LevelSetQuadrature.BernsteinPolynomial","page":"References","title":"LevelSetQuadrature.BernsteinPolynomial","text":"BernsteinPolynomial{D,T}\n\nTODO: docstring\n\n\n\n\n\n","category":"type"},{"location":"references/#LevelSetQuadrature.GradientDual","page":"References","title":"LevelSetQuadrature.GradientDual","text":"GradientDual{N,T}\n\nStructure used to reprensent the gradient of a function f : ℝᴺ → T for the purpose of performing forward-mode automatic differentiation.\n\n\n\n\n\n","category":"type"},{"location":"references/#LevelSetQuadrature.ImplicitDomain","page":"References","title":"LevelSetQuadrature.ImplicitDomain","text":"ImplicitDomain{N,T} <: AbstractDomain{N,T}\n\nDomain defined by\n\n\n\n\n\n","category":"type"},{"location":"references/#LevelSetQuadrature.LinearizationDual","page":"References","title":"LevelSetQuadrature.LinearizationDual","text":"struct LinearizationDual{D,T<:Real}\n\nLinearization of a D-dimensional function f : 𝐑ᴰ → T with a strict bound on the remainder.\n\nLinearizationDual objects are constructed given a function/functor f and a HyperRectangle rec using l = LinearizationDual(f,rec). The object l represents an approximation of f inside of rec in the following sense: |f(𝐱) - α - β⋅(𝐱 - 𝐱₀)| < ϵ ∀ 𝐱 ∈ rec, where l = LinearizationDual(f,rec) and α = value(l), β = gradient(l), ϵ = remainder(l).\n\n\n\n\n\n","category":"type"},{"location":"references/#LevelSetQuadrature.bound-Union{Tuple{D}, Tuple{Function, HyperRectangle{D}}} where D","page":"References","title":"LevelSetQuadrature.bound","text":"bound(f,rec::HyperRectangle)\n\nReturn an upper bound δ such that |f(x) - f(center(rec))| ≤ δ for x ∈ rec. By default, the upper bound is computed using a bounded LinearizationDual of f on rec.\n\nThis method can be overloaded for specific types typeof(f) if a more efficient way of computing the bound is known (e.g. if f is affine).\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.gradient-Tuple{LevelSetQuadrature.LinearizationDual}","page":"References","title":"LevelSetQuadrature.gradient","text":"gradient(l::Linearization)\n\nGradient of l at the center of the domain(l).\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.power2bernstein-Union{Tuple{Array{<:Real, D}}, Tuple{D}, Tuple{Array{<:Real, D}, HyperRectangle{D}}, Tuple{Array{<:Real, D}, HyperRectangle{D}, Any}} where D","page":"References","title":"LevelSetQuadrature.power2bernstein","text":"power2bernstein\n\nTODO: document and write a jldoctest\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.prune!-Tuple{LevelSetQuadrature.AbstractDomain}","page":"References","title":"LevelSetQuadrature.prune!","text":"prune!(Ω)\n\nPrune the functions specifying the domain Ω and return the CellType of the domain.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.quadgen-Union{Tuple{N}, Tuple{Function, HyperRectangle{N}, Any}, Tuple{Function, HyperRectangle{N}, Any, Any}} where N","page":"References","title":"LevelSetQuadrature.quadgen","text":"quadgen(ϕ,∇ϕ,U,s;order=5,maxdepth=20,maxgradient=20)\n\nTODO: extensively document this function\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.reference_cube-Union{Tuple{Val{D}}, Tuple{D}} where D","page":"References","title":"LevelSetQuadrature.reference_cube","text":"reference_cube(::Val{D})\n\nReturn the [0,1]ᴰ reference domain as a HyperRectangle{D,Float64}.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.remainder-Tuple{LevelSetQuadrature.LinearizationDual}","page":"References","title":"LevelSetQuadrature.remainder","text":"remainder(l::Linearization)\n\nAn upper bound for the remainder of l.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.sgn-NTuple{4, Any}","page":"References","title":"LevelSetQuadrature.sgn","text":"sgn(m,s,surface,side)\n\nCompute the sign of the upper and lower restrictions of a function ψ with sign s. The side variable is eithe 1 for the upper restriction, or -1 for the lower restriction, and m is the sign of ∇ψ ⋅ eₖ.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.tensorquad-Union{Tuple{D}, Tuple{HyperRectangle{D}, Any, Any}} where D","page":"References","title":"LevelSetQuadrature.tensorquad","text":"tensorquad(rec::HyperRectangle,x,w)\n\nGiven the one-dimensional quadrature nodes x and  weights w for integrating a function between 0 and 1, return the nodes and weights of the corresponding tensor-product quadrature on rec.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.value-Tuple{LevelSetQuadrature.LinearizationDual}","page":"References","title":"LevelSetQuadrature.value","text":"value(l::Linearization)\n\nValue of l at the center of the domain(l).\n\n\n\n\n\n","category":"method"},{"location":"references/#WavePropBase.domain-Tuple{LevelSetQuadrature.LinearizationDual}","page":"References","title":"WavePropBase.domain","text":"domain(l::Linearization) --> HyperRectangle\n\nDomain of validity of the linearization l.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Getting started","title":"Getting started","text":"CurrentModule = LevelSetQuadrature","category":"page"},{"location":"#home-section","page":"Getting started","title":"LevelSetQuadrature.jl","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"High-order quadratures for implicitly defined surfaces","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"note: Algorithm and references\nThis package implements, with minor modifications, the algorithms described in:Saye 2015\nSaye 2022If you find this work useful, please consider citing the original references.  TODO: properly cite.","category":"page"},{"location":"#Overview","page":"Getting started","title":"Overview","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"This package allows you to generate high-order quadrature nodes and weights for integrating regular functions over d-dimensional volumes/surfaces of the form:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Omega = boldsymbolx in U subset mathbbR^d quad st quad  phi(boldsymbolx)  0 quad textand quad Gamma = boldsymbolx in U subset mathbbR^d quad st quad  phi(boldsymbolx) = 0","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"where U is the bounding box described as a d-dimensional rectangle, and phi  mathbbR^d to mathbbR is a simple Julia function (more on this later). We refer to Omega and Gamma as implicit domains because they are defined implicitly through the function phi. For example, we may describe a disk through","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"using LevelSetQuadrature\nusing Plots\n\nf   = (x) -> x[1]^2 + x[2]^2 - 1\nlb, ub = -1.1, 1.1\n# plot\nxx = yy = lb:0.05:ub\nfig = contour(xx,yy,(x,y) -> f((x,y));levels=[0],aspect_ratio=:equal,cbar=nothing)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"To generate a quadrature for the disk we may then call quadgen as follows:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"rec = HyperRectangle((lb,lb),(ub,ub))\nx,w = quadgen(f,rec,:interior)\nprintln(\"error in area: \", abs(sum(w) - pi))","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"This produces the following set of nodes","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"scatter!(fig,PlotPoints(),x,label=\"interior nodes\",ms=3,m=:cross)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Similarly, a surface quadrature can be created by passing :surface instead of :interior as an argument to quadgen:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"x,w = quadgen(f,rec,:surface)\nprintln(\"error in perimeter:\", abs(sum(w) - 2*pi))","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"scatter!(fig,PlotPoints(),x,label=\"surface nodes\",ms=3)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"As long as the implicit function defining the domain is simple enough (e.g. a polynomial), changing the domain should be very straightforward:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"using LevelSetQuadrature\nusing Plots\n# Cassini oval: https://en.wikipedia.org/wiki/Cassini_oval\nconst a = 1.1\nconst c = 1.0\nf = (x) -> (x[1]^2 + x[2]^2)^2 - 2*c^2*(x[1]^2 - x[2]^2) - (a^4 - c^4)\n# a bounding box\nl,u = -1.5,1.5\nrec = HyperRectangle((l,l),(u,u))\n# generate the quadrature\nx,w = quadgen(f,rec,:interior)\n\n# plot the domain and quadrature\nxx = yy = l:0.05:u\ncontour(xx,yy,(x,y) -> f((x,y));levels=[0],aspect_ratio=:equal,cbar=nothing,lw=3)\nscatter!(PlotPoints(),x,label=nothing,ms=3)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"There is nothing particularly important about dimension 2 (besides the convenient fact that 2d is easy to visualize) – the algorithms implemented are dimension agnostic. Computing a surface quadrature for the unit sphere for instance can be achieved through:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"using LevelSetQuadrature\nusing Plots\nf = (x) -> sum(x.*x) - 1\nl,u = -1.5,1.5\nrec = HyperRectangle((l,l,l),(u,u,u))\nx,w = quadgen(f,rec,:surface)\nprintln(\"error in surface area: \", sum(w) - 4π)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"# plot the nodes\nscatter(PlotPoints(),x,label=nothing,ms=3,m=:cross)","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"There are several keyword arguments that can be passed to quadgen which provide further control on the generated quadrature; you should consult the quadgen documentation for more details.","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"warning: Limitation to simple functions\nIn order to build a high-order quadrature, the algorithm implemented in this package must build rigorous bounds on the implicit function passed and its derivatives. For that purpose, we use the bounded linearization idea of Saye, implemented in the LinearizationDual type, which performs an operation akin to automatic differentiation in order to propagate the bounds on the linearization of functions. Because we have only implemented a small subset of the \"linearization rules\", it is possible/likely that the call quadgen(f,rec) will fail if the algorithm does not manage to propagate the bounds on f and its gradient; that is, if the call f(::LinearizationDual) fails. Feel free to open an issue (or a PR) if quadgen does not work with your simple function.","category":"page"},{"location":"#Polynomial-domains","page":"Getting started","title":"Polynomial domains","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"In the case of polynomial domains ...","category":"page"},{"location":"#Advanced-usage","page":"Getting started","title":"Advanced usage","text":"","category":"section"},{"location":"#Other-limitations-and-known-issues","page":"Getting started","title":"Other limitations and known issues","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"Here are some know limitations of the package:","category":"page"},{"location":"","page":"Getting started","title":"Getting started","text":"Automatic differentiation\nSingular points\n","category":"page"},{"location":"#Index","page":"Getting started","title":"Index","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"","category":"page"}]
}
