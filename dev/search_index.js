var documenterSearchIndex = {"docs":
[{"location":"references/#references-section","page":"References","title":"References","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Modules = [LevelSetQuadrature]","category":"page"},{"location":"references/#LevelSetQuadrature.Linearization","page":"References","title":"LevelSetQuadrature.Linearization","text":"struct Linearization{D,T<:Real}\n\nLinearization of a D-dimensional function f : 𝐑ᴰ → T with a strict bound on the remainder.\n\nLinearization objects are constructed given a function/functor f and a HyperRectangle rec using l = Linearization(f,rec). The object l represents an approximation of f inside of rec in the following sense: |f(𝐱) - α - β⋅(𝐱 - 𝐱₀)| < ϵ ∀ 𝐱 ∈ rec, where l = Linearization(f,rec) and α = value(l), β = gradient(l), ϵ = remainder(l)\n\n\n\n\n\n","category":"type"},{"location":"references/#LevelSetQuadrature.bound-Union{Tuple{D}, Tuple{Function, HyperRectangle{D}}} where D","page":"References","title":"LevelSetQuadrature.bound","text":"bound(f,rec::HyperRectangle)\n\nReturn an upper bound δ such that |f(x) - f(center(rec))| ≤ δ for x ∈ rec. By default, the upper bound is computed using a bounded Linearization of f on rec.\n\nThis method can be overloaded for specific types typeof(f) if a more efficient way of computing the bound is known (e.g. if f is affine).\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.gradient-Tuple{Linearization}","page":"References","title":"LevelSetQuadrature.gradient","text":"value(l::Linearization)\n\nGradient of l at the center of the domain(l).\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.remainder-Tuple{Linearization}","page":"References","title":"LevelSetQuadrature.remainder","text":"remainder(l::Linearization)\n\nAn upper bound for the remainder of l.\n\n\n\n\n\n","category":"method"},{"location":"references/#LevelSetQuadrature.value-Tuple{Linearization}","page":"References","title":"LevelSetQuadrature.value","text":"value(l::Linearization)\n\nValue of l at the center of the domain(l).\n\n\n\n\n\n","category":"method"},{"location":"references/#WavePropBase.domain-Tuple{Linearization}","page":"References","title":"WavePropBase.domain","text":"domain(l::Linearization) --> HyperRectangle\n\nDomain of validity of the linearization l.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Getting started","title":"Getting started","text":"CurrentModule = LevelSetQuadrature","category":"page"},{"location":"#home-section","page":"Getting started","title":"LevelSetQuadrature.jl","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"High-order quadratures for implicitly defined surfaces","category":"page"},{"location":"#Overview","page":"Getting started","title":"Overview","text":"","category":"section"},{"location":"#Index","page":"Getting started","title":"Index","text":"","category":"section"},{"location":"","page":"Getting started","title":"Getting started","text":"","category":"page"}]
}
