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

function quadgen(p::DynamicPolynomials.Polynomial,U,s::Symbol;kwargs...)
    signs,surf = _symbol_to_signs_and_surf(s)
    coefs = monomial_coefs(p)
    ϕ     = power2bernstein(coefs,U)
    Ω     = BernsteinDomain([ϕ],signs)
    quadgen(Ω,surf;kwargs...)
end
