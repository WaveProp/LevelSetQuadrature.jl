struct Linearization{N,T}
    α::T
    β::SVector{N,T}
    δ::SVector{N,T}
    ϵ::T
end

function Base.:+(u::Linearization,v::Linearization)
    @assert u.δ == v.δ
    Linearization(u.α + v.α, u.β + v.β, u.δ, u.ϵ + v.ϵ)
end

function Base.:-(u::Linearization,v::Linearization)
    @assert u.δ == v.δ
    Linearization(u.α - v.α, u.β - v.β, u.δ, u.ϵ + v.ϵ)
end

function Base.:*(u::Linearization,v::Linearization)
    @assert u.δ == v.δ
    # FIXME
end
