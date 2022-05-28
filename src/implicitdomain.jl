struct ImplicitDomain
    Ψ::Vector{Function}
    ∇Ψ::Vector{Function}
    signs::Vector{Int}
    rec::HyperRectangle
end
