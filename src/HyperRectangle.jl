struct HyperRectangle{D,T<:Real}
    L::SVector{D,T}
    U::SVector{D,T}
    level::Int64
end

HyperRectangle(L, U) = HyperRectangle(L, U, 0)

center(rec::HyperRectangle) = (rec.U + rec.L) / 2
δ(rec::HyperRectangle) = (rec.U - rec.L) / 2

function divide(rec::HyperRectangle{D}) where{D}
    ax = argmax(rec.U .- rec.L)
    U1 = Vector(rec.U); U1[ax] = (rec.L[ax] + rec.U[ax])/2
    L2 = Vector(rec.L); L2[ax] = (rec.L[ax] + rec.U[ax])/2
    HyperRectangle(rec.L, SVector{D}(U1), rec.level), HyperRectangle(SVector{D}(L2), rec.U, rec.level)
end

function section(rec::HyperRectangle{D}, ax::Integer) where{D}
    @assert D ≥ ax
    HyperRectangle(deleteat(rec.L, ax), deleteat(rec.U, ax), rec.level)
end
