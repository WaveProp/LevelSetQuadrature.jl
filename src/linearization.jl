struct Linearization{D,T<:Real}
    rec::HyperRectangle{D,T}
    α::T
    β::SVector{D,T}
    δ::SVector{D,T}
    ϵ::T
end

"""
    value(l::Linearization)

Return the value of `l` at the center of `rectangle(l)`.
"""
value(l::Linearization)  = l.α

gradient(l::Linearization)  = l.β

remainder(l::Linearization)  = l.ϵ

domain(l::Linearization) = l.rec

half_width(l::Linearization) = l.δ

"""
    bound(l::Linearization)

Maximum difference between ...
"""
function bound(l::Linearization)
    δ = half_width(l)
    β = gradient(l)
    dot(abs.(β),δ) + remainder(l)
end

## Constructors
Linearization(rec::HyperRectangle, α::T, β::SVector{D,T}, ϵ::T) where {D,T} = Linearization{D,T}(rec, α, β, half_width(rec), ϵ)

## addition
Base.:+(u::Linearization{D,T}, c::Real) where{D,T} = Linearization{D,T}(u.rec, u.α+c, u.β, u.δ, u.ϵ)
Base.:+(c::Real, u::Linearization{D,T}) where{D,T} = Linearization{D,T}(u.rec, u.α+c, u.β, u.δ, u.ϵ)
function Base.:+(u::Linearization{D,T},v::Linearization{D,T}) where{D,T}
    # @assert u.rec == v.rec
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
    # @assert u.rec == v.rec
    l1 = abs.(u.β)' * u.δ
    l2 = abs.(v.β)' * u.δ
    Linearization{D,T}(u.rec, u.α * v.α, u.α*v.β + v.α*u.β, u.δ, l1*l2 + (abs(u.α)+l1)*v.ϵ + (abs(v.α)+l2)*u.ϵ + u.ϵ*v.ϵ)
end

## division
Base.:/(u::Linearization{D,T}, c::Real) where{D,T} = Linearization{D,T}(u.rec, u.α/c, u.β./c, u.δ, u.ϵ/abs(c))

# TODO: explain this function
function dual_variables(rec::HyperRectangle{D,T}) where {D,T}
    xc  = center(rec)
    x̂ = ntuple(D) do dim
        β = ntuple(i->i==dim ? one(T) : zero(T),D) |> SVector
        Linearization(rec,xc[dim],β,half_width(rec),zero(T))
    end
    return SVector{D, Union{Linearization{D,T}, T}}(x̂)
end

function linearization(f,rec::HyperRectangle{D}) where {D}
    x̂ = dual_variables(rec)
    f(x̂)
end

function Base.:^(l::Linearization,p::Integer)
    @assert p ≥ 1
    if p == 1
        return l
    else
        l*(l^(p-1))
    end
end

function bound(f::Function, rec::HyperRectangle)
    # return 0
    f̂ = linearization(f,rec)
    bound(f̂)
end

function ∇(f::Function, D, h=1e-5)
    grad(x) = ntuple(D) do d
        x̃ = [i == d ? x[i]+h : x[i] for i = 1:D]
        (f(x̃) - f(x)) / h
    end |> SVector
end


# function Base.show(io::IO,l::Linearization)
#     println("Linearization over ($rec) with $(value(l)) + $(gradient(l))⋅(x-xc) + [-$(remainder(l)),$(remainder(l))]")
# end

# finite difference gradient implementation.
function gradient(f::Function,x::SVector{N,T}) where {N,T}
    h = (eps(T))^(1/3)
    ntuple(N) do d
        xp = ntuple(i-> i==d ? x[i]+h : x[i] , N) |> SVector
        xm = ntuple(i-> i==d ? x[i]-h : x[i] , N) |> SVector
        (f(xp) - f(xm))/(2h)
    end |> SVector
end

Base.eps(::Type{Linearization{N,T}}) where {N,T} = eps(T)
