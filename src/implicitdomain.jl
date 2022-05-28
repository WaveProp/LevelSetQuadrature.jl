@enum CellType empty_cell whole_cell cut_cell

struct ImplicitDomain{N,T}
    Ψ::Vector{<:Function}
    ∇Ψ::Vector{<:Function}
    signs::Vector{Int}
    rec::HyperRectangle{N,T}
    celltype::CellType
    function ImplicitDomain(Ψ,∇Ψ,signs,rec::HyperRectangle{N,T}) where {N,T}
        ctype = _prune!(Ψ,∇Ψ,rec,signs)
        new{N,T}(Ψ,∇Ψ,signs,rec::HyperRectangle{N,T},ctype)
    end
end

cell_type(Ω::ImplicitDomain) = Ω.celltype

"""
    prune!(Ω)

Prune the functions specifying the domain `Ω` and return the `CellType` of the domain.
"""
function prune!(Ω::ImplicitDomain)
    _prune!(Ω.Ψ,Ω.∇Ψ,Ω.rec,Ω.signs)
end

function _prune!(Ψ,∇Ψ,rec,signs)
    delInd = Vector{Int}()
    for (i, ψ) in enumerate(Ψ)
        si = signs[i]
        t = cell_type(ψ,si,rec)
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

function cell_type(ψ,s,rec)
    δ = bound(ψ, rec)
    ψc = ψ(center(rec))
    abs(ψc) ≥ δ || (return cut_cell)
    if s * ψc ≥ 0
        # intersection is the whole rec
        return whole_cell
    else
        # intersection is empty, return immediately
        return empty_cell
    end
end
