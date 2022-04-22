struct HyperRectangle{D,T<:Real}
    L::SVector{D,T}
    U::SVector{D,T}
    level::Int64
end

HyperRectangle(L, U) = HyperRectangle(L, U, 0) # build root

center(rec::HyperRectangle)     = (rec.U + rec.L) / 2
width(rec::HyperRectangle)      = (rec.U - rec.L)
half_width(rec::HyperRectangle) = width(rec) / 2


function divide(rec::HyperRectangle{D}) where{D}
    ax = argmax(rec.U .- rec.L)
    m = (rec.L[ax] + rec.U[ax])/2
    U1 = ntuple(dim -> dim == ax ? m : rec.U[dim], D) |> SVector
    L2 = ntuple(D) do dim
        dim == ax ? m : rec.L[dim]
    end |> SVector
    HyperRectangle(rec.L, U1, rec.level+1), HyperRectangle(L2, rec.U, rec.level+1)
end

function section(rec::HyperRectangle{D}, ax::Integer) where{D}
    @assert 1 ≤ ax ≤ D
    HyperRectangle(deleteat(rec.L, ax), deleteat(rec.U, ax), rec.level)
end

function Base.show(rec::HyperRectangle{D}) where{D}
    str = "[$(rec.L[1]), $(rec.U[1])]"
    for d = 2:D 
        str *= "×[$(rec.L[d]), $(rec.U[d])]"
    end
    println(str)
end