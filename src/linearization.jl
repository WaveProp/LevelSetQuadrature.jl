struct Linearization{D,T<:Real}
    rec::HyperRectangle{D,T}
    α::T
    β::SVector{D,T}
    δ::SVector{D,T}
    ϵ::T
end

## Constructors
Linearization(rec::HyperRectangle, α::T, β::SVector{D,T}, ϵ::T) where{D,T} = Linearization{D,T}(rec, α, β, δ(rec), ϵ)

## addition
Base.:+(u::Linearization{D,T}, c::Real) where{D,T} = Linearization{D,T}(u.rec, u.α+c, u.β, u.δ, u.ϵ)
Base.:+(c::Real, u::Linearization{D,T}) where{D,T} = Linearization{D,T}(u.rec, u.α+c, u.β, u.δ, u.ϵ)
function Base.:+(u::Linearization{D,T},v::Linearization{D,T}) where{D,T}
    @assert u.rec == v.rec
    Linearization{D,T}(u.rec, u.α + v.α, u.β + v.β, u.δ, u.ϵ + v.ϵ)
end

## subtraction
Base.:-(u::Linearization{D,T}) where{D,T} = Linearization{D,T}(u.rec, -u.α, -u.β, u.δ, u.ϵ)
Base.:-(u::Linearization{D,T}, c::Real) where{D,T} = Linearization{D,T}(u.rec, u.α-c, u.β, u.δ, u.ϵ)
Base.:-(c::Real, u::Linearization{D,T}) where{D,T} = Linearization{D,T}(u.rec, c-u.α, -u.β, u.δ, u.ϵ)
function Base.:-(u::Linearization{D,T},v::Linearization{D,T}) where{D,T}
    @assert u.rec == v.rec
    Linearization{D,T}(u.rec, u.α - v.α, u.β - v.β, u.δ, u.ϵ + v.ϵ)
end

## multiplication
Base.:*(u::Linearization{D,T}, c::Real) where{D,T} = Linearization{D,T}(u.rec, u.α*c, u.β*c, u.δ, u.ϵ*abs(c))
Base.:*(c::Real, u::Linearization{D,T}) where{D,T} = Linearization{D,T}(u.rec, u.α*c, u.β*c, u.δ, u.ϵ*abs(c))
function Base.:*(u::Linearization{D,T},v::Linearization{D,T}) where{D,T}
    @assert u.rec == v.rec
    l1 = abs.(u.β)' * u.δ
    l2 = abs.(v.β)' * u.δ
    Linearization{D,T}(u.rec, u.α * v.α, u.α*v.β + v.α*u.β, u.δ, l1*l2 + (abs(u.α)+l1)*v.ϵ + (abs(v.α)+l2)*u.ϵ + u.ϵ*v.ϵ)
end

function bound(f::Function, U::HyperRectangle)
    0# replace by max/min (extrema) of f on a grid of U
    # FIXME 
end