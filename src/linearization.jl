abstract type AbstractDual{N,T} end

"""
    GradientDual{N,T}

Structure used to reprensent the gradient of a function `f : ℝᴺ → T` for the
purpose of performing forward-mode automatic differentiation.
"""
struct GradientDual{N,T} <: AbstractDual{N,T}
    α::T
    β::SVector{N,T}
end

value(l::GradientDual)    = l.α
gradient(l::GradientDual) = l.β

"""
    struct LinearizationDual{D,T<:Real}

Linearization of a `D`-dimensional function `f : 𝐑ᴰ → T` with a strict bound on
the remainder.

`LinearizationDual` objects are constructed given a function/functor `f` and a
`HyperRectangle` `rec` using `l = LinearizationDual(f,rec)`. The object `l`
represents an approximation of `f` inside of `rec` in the following sense:
`|f(𝐱) - α - β⋅(𝐱 - 𝐱₀)| < ϵ ∀ 𝐱 ∈ rec`, where `l = LinearizationDual(f,rec)` and `α =
value(l)`, `β = gradient(l)`, `ϵ = remainder(l)`.
"""
struct LinearizationDual{D,T} <: AbstractDual{D,T}
    α::T
    β::SVector{D,T}
    ϵ::T
    rec::HyperRectangle{D,T}
end

"""
    value(l::Linearization)

Value of `l` at the center of the `domain(l)`.
"""
value(l::LinearizationDual) = l.α

"""
    gradient(l::Linearization)

Gradient of `l` at the center of the `domain(l)`.
"""
gradient(l::LinearizationDual) = l.β

"""
    remainder(l::Linearization)

An upper bound for the remainder of `l`.
"""
remainder(l::LinearizationDual) = l.ϵ

"""
    domain(l::Linearization) --> HyperRectangle

Domain of validity of the linearization `l`.
"""
domain(l::LinearizationDual) = l.rec

## addition
Base.:+(u::GradientDual, c::Real) = GradientDual(value(u) + c, gradient(u))
Base.:+(c::Real, u::GradientDual) = u + c
Base.:+(u::LinearizationDual, c::Real) = LinearizationDual(value(u) + c, gradient(u), remainder(u), domain(u))
Base.:+(c::Real, u::LinearizationDual) = u + c
function Base.:+(u::LinearizationDual{N,T}, v::LinearizationDual{N,T}) where {N,T}
    rec = _result_domain(u,v)
    LinearizationDual(value(u) + value(v), gradient(u) + gradient(v), u.ϵ + v.ϵ, rec)
end
function Base.:+(u::GradientDual, v::GradientDual)
    GradientDual(value(u) + value(v), gradient(u) + gradient(v))
end

## subtraction
Base.:-(u::LinearizationDual) = (-1) * u
Base.:-(u::GradientDual) = (-1) * u
Base.:-(u::LinearizationDual, c::Real) = u + (-c)
Base.:-(u::GradientDual, c::Real) = u + (-c)
Base.:-(c::Real, u::LinearizationDual) = c + (-u)
Base.:-(c::Real, u::GradientDual) = c + (-u)
Base.:-(u::LinearizationDual, v::LinearizationDual) = u + (-v)
Base.:-(u::GradientDual, v::GradientDual) = u + (-v)

## multiplication
Base.:*(u::LinearizationDual, c::Real) = LinearizationDual(value(u) * c, gradient(u) * c, abs(c) * remainder(u), domain(u))
Base.:*(u::GradientDual, c::Real) = GradientDual(value(u) * c, gradient(u) * c)
Base.:*(c::Real, u::LinearizationDual) = u * c
Base.:*(c::Real, u::GradientDual) = u * c
function Base.:*(u::LinearizationDual, v::LinearizationDual)
    rec = _result_domain(u,v)
    δ = half_width(u)
    l1 = dot(abs.(gradient(u)), δ)
    l2 = dot(abs.(gradient(v)), δ)
    LinearizationDual(value(u) * value(v), value(u) * gradient(v) + value(v) * gradient(u), l1 * l2 + (abs(value(u)) + l1) * remainder(v) + (abs(value(v)) + l2) * remainder(u) + remainder(u) * remainder(v), rec)
end

# since `LinearizationDual` behaves like a number, multiplying a `LinearizationDual` by
# a vector of the same type is handled like scalar*vector
function Base.:*(α::LinearizationDual{N, T}, β::SVector{M, LinearizationDual{N, T}}) where {N,M,T}
    map(b -> α*b,β)
end

function Base.:*(u::GradientDual, v::GradientDual)
    GradientDual(value(u) * value(v), value(u) * gradient(v) + value(v) * gradient(u))
end

# division
Base.:/(u::LinearizationDual, c::Real) = LinearizationDual(value(u) / c, gradient(u) / c, remainder(u) / abs(c), domain(u))
Base.:/(u::GradientDual, c::Real) = GradientDual(value(u) / c, gradient(u) / c)

# dual basis
function linearization_basis(rec::HyperRectangle{D,T}) where {D,T}
    xc = center(rec)
    x̂ = svector(D) do dim
        β = svector(i -> i == dim ? one(T) : zero(T), D)
        LinearizationDual(xc[dim], β, zero(T), rec)
    end
    return x̂
end

function gradient_basis(x::SVector{D,T}) where {D,T}
    x̂ = svector(D) do dim
        β = svector(i -> i == dim ? one(T) : zero(T), D)
        GradientDual(x[dim], β)
    end
    return x̂
end

# TODO: document
function linearization(f, rec::HyperRectangle)
    x̂ = linearization_basis(rec)
    f(x̂)
end

# TODO: document
function gradient(f, x::SVector)
    x̂ = gradient_basis(x)
    gradient(f(x̂))
end

# power
function Base.:^(l::LinearizationDual, p::Integer)
    @assert p ≥ 1
    if p == 1
        return l
    else
        l * (l^(p - 1))
    end
end

function Base.:^(l::GradientDual, p::Integer)
    @assert p ≥ 1
    if p == 1
        return l
    else
        l * (l^(p - 1))
    end
end

# TODO: add cos/sin

"""
    bound(f,rec::HyperRectangle)

Return an upper bound `δ` such that `|f(x) - f(center(rec))| ≤ δ` for `x ∈ rec`.
By default, the upper bound is computed using a bounded [`LinearizationDual`](@ref) of `f`
on `rec`.

This method can be overloaded for specific types `typeof(f)` if a more efficient
way of computing the bound is known (e.g. if `f` is affine).
"""
function bound(f::Function, rec::HyperRectangle{D}) where {D}
    f̂ = linearization(f, rec)
    bound(f̂)
end

#
function approximate_bound(f::Function, rec::HyperRectangle{D}, n=10) where {D}
    sz = ntuple(i -> n, D)
    msh = UniformCartesianMesh(rec, sz)
    iter = NodeIterator(msh)
    fmin, fmax = extrema(f, iter)
    return fmin,fmax
end

# HACK: automatic definition of gradient as required by this package (i.e. as an SVector
# of functions). Maybe not the most efficient one, but does not seem to affect
# the performance in any significant way. A better way could be to get the
# partial independently
function gradient(ϕ::Function,::Val{D}) where {D}
    f = (x) -> gradient(ϕ,x)
    svector(i-> x -> f(x)[i],D)
end

half_width(l::LinearizationDual) = half_width(domain(l))

function bound(l::LinearizationDual)
    α = value(l)
    δ = half_width(l)
    β = gradient(l)
    Δ = dot(abs.(β), δ) + remainder(l)
    α - Δ, α + Δ
end

bound(l::SVector{<:Any,<:LinearizationDual}) = bound.(l) # for gradients

# required by `gradient_basis` to autodiff through `LinearizationDual`
function Base.one(::Type{LinearizationDual{N, T}}) where {N,T}
    α = one(T)
    β = svector(i->zero(T),N)
    ϵ = zero(T)
    rec = HyperRectangle(β,β)
    LinearizationDual(α,β,ϵ,rec)
end
function Base.zero(::Type{LinearizationDual{N, T}}) where {N,T}
    α = zero(T)
    β = svector(i->zero(T),N)
    ϵ = zero(T)
    rec = HyperRectangle(β,β)
    LinearizationDual(α,β,ϵ,rec)
end

# domain of operation between two linearizations
function _result_domain(u::LinearizationDual{N,T},v::LinearizationDual{N,T}) where {N,T}
    # HACK: the zero rec is used when building one(::Type{LinearizationDual})
    # and zero(::Type{LinearizationDual}), so we must check those before
    # asserting equality
    rec = HyperRectangle(svector(i->zero(T),N),svector(i->zero(T),N))
    if domain(u) == rec
        return domain(v)
    elseif domain(v) == rec
        return domain(u)
    else
        @assert domain(u) == domain(v) "domains must be identical"
        return domain(u)
    end
end
