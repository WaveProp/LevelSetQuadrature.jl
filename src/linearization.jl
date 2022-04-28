"""
    struct Linearization{D,T<:Real}

Linearization of a `D`-dimensional function `f : 𝐑ᴰ → T` with a strict bound on
the remainder.

`Linearization` objects are constructed given a function/functor `f` and a
`HyperRectangle` `rec` using `l = Linearization(f,rec)`. The object `l`
represents an approximation of `f` inside of `rec` in the following sense:
`|f(𝐱) - α - β⋅(𝐱 - 𝐱₀)| < ϵ ∀ 𝐱 ∈ rec`, where `l = Linearization(f,rec)` and `α =
value(l)`, `β = gradient(l)`, `ϵ = remainder(l)`
"""
struct Linearization{D,T<:Real}
    α::T
    β::SVector{D,T}
    ϵ::T
    rec::HyperRectangle{D,T}
end

"""
    value(l::Linearization)

Value of `l` at the center of the `domain(l)`.
"""
value(l::Linearization) = l.α

"""
    value(l::Linearization)

Gradient of `l` at the center of the `domain(l)`.
"""
gradient(l::Linearization) = l.β

"""
    remainder(l::Linearization)

An upper bound for the remainder of `l`.
"""
remainder(l::Linearization) = l.ϵ

"""
    domain(l::Linearization) --> HyperRectangle

Domain of validity of the linearization `l`.
"""
WavePropBase.domain(l::Linearization) = l.rec

WavePropBase.half_width(l::Linearization) = half_width(domain(l))

function bound(l::Linearization)
    δ = half_width(l)
    β = gradient(l)
    dot(abs.(β), δ) + remainder(l)
end

## addition
Base.:+(u::Linearization, c::Real) = Linearization(value(u) + c, gradient(u), remainder(u), domain(u))

Base.:+(c::Real, u::Linearization) = u + c

function Base.:+(u::Linearization, v::Linearization)
    @assert domain(u) == domain(v) "domains must be identical"
    Linearization(value(u) + value(v), gradient(u) + gradient(v), u.ϵ + v.ϵ, domain(u))
end

## subtraction

Base.:-(u::Linearization) = (-1) * u

Base.:-(u::Linearization, c::Real) = u + (-c)

Base.:-(c::Real, u::Linearization) = c + (-u)

Base.:-(u::Linearization, v::Linearization) = u + (-v)

## multiplication
Base.:*(u::Linearization, c::Real) = Linearization(value(u) * c, gradient(u) * c, abs(c) * remainder(u), domain(u))

Base.:*(c::Real, u::Linearization) = u * c

function Base.:*(u::Linearization, v::Linearization)
    @assert domain(u) == domain(v)
    δ = half_width(u)
    l1 = dot(abs.(gradient(u)), δ)
    l2 = dot(abs.(gradient(v)), δ)
    Linearization(value(u) * value(v), value(u) * gradient(v) + value(v) * gradient(u), l1 * l2 + (abs(value(u)) + l1) * remainder(v) + (abs(value(v)) + l2) * remainder(u) + remainder(u) * remainder(v), domain(u))
end

## division
Base.:/(u::Linearization, c::Real) = Linearization(value(u) / c, gradient(u) / c, remainder(u) / abs(c), domain(u))

# TODO: explain this function
function dual_variables(rec::HyperRectangle{D,T}) where {D,T}
    xc = center(rec)
    x̂ = ntuple(D) do dim
        β = svector(i -> i == dim ? one(T) : zero(T), D)
        Linearization(xc[dim], β, zero(T), rec)
    end
    # FIXME: why this? I suspect it is related to the hack in the `V` function
    return SVector{D, Union{Linearization{D,T}, T}}(x̂)
end

function linearization(f, rec::HyperRectangle)
    x̂ = dual_variables(rec)
    f(x̂)
end

function Base.:^(l::Linearization, p::Integer)
    @assert p ≥ 1
    if p == 1
        return l
    else
        l * (l^(p - 1))
    end
end

"""
    bound(f,rec::HyperRectangle)

Return an upper bound `δ` such that `|f(x) - f(center(rec))| ≤ δ` for `x ∈ rec`.
By default, the upper bound is computed using a bounded [`Linearization`](@ref) of `f`
on `rec`.

This method can be overloaded for specific types `typeof(f)` if a more efficient
way of computing the bound is known (e.g. if `f` is affine).
"""
function bound(f::Function, rec::HyperRectangle{D}) where {D}
    # approximate_bound(f,rec)
    f̂ = linearization(f, rec)
    bound(f̂)
end

#
function approximate_bound(f::Function, rec::HyperRectangle{D}, n=100) where {D}
    sz = ntuple(i -> n, D)
    msh = UniformCartesianMesh(rec, sz)
    iter = NodeIterator(msh)
    fmin, fmax = extrema(f, iter)
    fc = f(center(rec))
    return max(fmax - fc, fc - fmin)
end
