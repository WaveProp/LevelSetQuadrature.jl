#=
    File including wrappers to handle polynomials from `DynamicPolynomials`.
    This is conditionally loaded if `DynamicPolynomials` is loaded before `using
    LevelSetQuadrature`.
=#

import .DynamicPolynomials

function monomial_coefs(p::DynamicPolynomials.Polynomial)
    T = DynamicPolynomials.coefficienttype(p)
    sz    = mapreduce(m->DynamicPolynomials.exponents(m),(a,b)->max.(a,b),DynamicPolynomials.monomials(p)) .+ 1
    coefs = zeros(T,sz...)
    for t in p
        m = DynamicPolynomials.exponents(t)
        v = DynamicPolynomials.coefficient(t)
        I = m .+ 1
        coefs[Tuple(I)...] = v
    end
    coefs
end

function quadgen!(X, W, p::DynamicPolynomials.Polynomial,U::HyperRectangle{N},s::Symbol;kwargs...) where {N}
    coefs = monomial_coefs(p)
    ϕ     = power2bernstein(coefs, U)
    quadgen!(X, W, ϕ, s;kwargs...)
end

function quadgen(p::DynamicPolynomials.Polynomial,U::HyperRectangle{N},s::Symbol;kwargs...) where {N}
    X = Vector{SVector{N,Float64}}()
    W = Vector{Float64}()
    quadgen!(X, W, p, U, s;kwargs...)
end
