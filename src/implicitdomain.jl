@enum CellType empty_cell whole_cell cut_cell

struct ImplicitDomain{N,T}
    Ψ::Vector{<:Function}
    ∇Ψ::Vector{SVector{N,<:Function}}
    signs::Vector{Int}
    rec::HyperRectangle{N,T}
    celltype::CellType
    function ImplicitDomain(Ψ, ∇Ψ, signs, rec::HyperRectangle{N,T}) where {N,T}
        ctype = _prune!(Ψ, ∇Ψ, rec, signs)
        new{N,T}(Ψ, ∇Ψ, signs, rec::HyperRectangle{N,T}, ctype)
    end
end

cell_type(Ω::ImplicitDomain) = Ω.celltype

"""
    prune!(Ω)

Prune the functions specifying the domain `Ω` and return the `CellType` of the domain.
"""
function prune!(Ω::ImplicitDomain)
    _prune!(Ω.Ψ, Ω.∇Ψ, Ω.rec, Ω.signs)
end

function _prune!(Ψ, ∇Ψ, rec, signs)
    delInd = Vector{Int}()
    for (i, ψ) in enumerate(Ψ)
        si = signs[i]
        t = cell_type(ψ, si, rec)
        if t == whole_cell
            # intersection is the whole rec, so ψ can be prune
            # @info "Whole cell"
            append!(delInd, i)
        elseif t == empty_cell
            # intersection is empty, return immediately
            return empty_cell
        end
    end
    deleteat!(signs, delInd)
    deleteat!(Ψ, delInd)
    deleteat!(∇Ψ, delInd)
    isempty(Ψ) && return whole_cell
    return cut_cell
end

function cell_type(ψ, s, rec)
    l, u = bound(ψ, rec)
    ψc = ψ(center(rec))
    l * u ≥ 0 || (return cut_cell)
    if s * ψc ≥ 0
        # intersection is the whole rec
        return whole_cell
    else
        # intersection is empty, return immediately
        return empty_cell
    end
end

function restrict(Ω::ImplicitDomain{N,T}, k, surf) where {N,T}
    Ψ = Ω.Ψ
    ∇Ψ = Ω.∇Ψ
    signs = Ω.signs
    rec = Ω.rec
    xc = center(rec)
    Ψ̃  = Function[]
    ∇Ψ̃ = SVector{N-1,Function}[]
    new_signs = empty(signs)
    for (ψ, s, ∇ψ) in zip(Ψ, signs, ∇Ψ)
        # why bound? dont we know that ∇Ψ[k] has a fixed sign on direction k?
        # pos_neg = bound(∇ψ[k],rec)[1] > 0 ? 1 : -1
        pos_neg = ∇ψ[k](xc) > 0 ? 1 : -1 # use sign?
        ψL = lower_restrict(ψ, rec, k)
        sL = sgn(pos_neg, s, surf, -1)
        ∇ψL = lower_restrict_grad(∇ψ, rec, k)
        ψU = upper_restrict(ψ, rec, k)
        sU = sgn(pos_neg, s, surf, 1)
        ∇ψU = upper_restrict_grad(∇ψ, rec, k)
        append!(Ψ̃, (ψL, ψU))
        append!(new_signs, (sL, sU))
        append!(∇Ψ̃, (∇ψL, ∇ψU))
    end
    Ω̃ = ImplicitDomain(Ψ̃, ∇Ψ̃, new_signs, section(rec, k))
    return Ω̃
end


function lower_restrict(ψ::Function,rec,k)
    a = low_corner(rec)[k]
    x -> ψ(insert(x, k, a))
end

function upper_restrict(ψ::Function,rec,k)
    a = high_corner(rec)[k]
    x -> ψ(insert(x, k, a))
end

function lower_restrict_grad(∇ψ::SVector{N}, rec, k) where {N}
    a = low_corner(rec)[k]
    ∇ψ′ = deleteat(∇ψ,k)
    (x) -> ∇ψ′(insert(x, k, a))
    svector(d -> (x) -> ∇ψ′[d](insert(x, k, a)), N-1)
end

function upper_restrict_grad(∇ψ::SVector{N}, rec, k) where {N}
    a   = high_corner(rec)[k]
    ∇ψ′ = deleteat(∇ψ,k)
    svector(d -> (x) -> ∇ψ′[d](insert(x, k, a)), N-1)
end
