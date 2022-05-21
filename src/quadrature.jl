function StaticArrays.insert(x̂::SVector{<:Any,<:Linearization{N,T}},k,y::T) where {N,T}
    rec = first(x̂) |> domain
    ŷ = Linearization(y,zero(SVector{N,T}),zero(T),rec)
    insert(x̂,k,ŷ)
end

function StaticArrays.insert(x::Real,k,y::Real)
    insert(SVector(x),k,y)
end

# TODO: @dongchen explain this function?
function sgn(m, s, S, σ)
    if m * s * σ > 0 || S
        if m * σ > 0
            return 1
        else
            return -1
        end
    else
        return 0
    end
end

"""
    tensorquad(rec,x,w)

Given a one-dimensional
"""
function tensorquad(rec::HyperRectangle{D}, GX, GW) where {D}
    xl, xu = low_corner(rec), high_corner(rec)
    if D == 1
        nodes = (GX .+ 1.0) ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1]) .+ low_corner(rec)[1]
        weights = GW ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1])
        return nodes, weights
    else
        rec1 = section(rec, 1)
        nodes = Vector{SVector{D,Float64}}()
        weights = Vector{Float64}()
        X, W = tensorquad(rec1, GX, GW)
        for (x, w) in zip(X, W)
            Y = (GX .+ 1.0) ./ 2.0 * (high_corner(rec)[1] - low_corner(rec)[1]) .+ low_corner(rec)[1]
            Ω = GW ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1])
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, 1, y))
                push!(weights, w * ω)
            end
        end
    end
    nodes, weights
end

function dim1NodesWeights(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, L, U, q, multi_zeros=true)
    roots = [L, U]
    if multi_zeros
        for ψ in Ψ
            union!(roots, find_zeros(ψ, L, U))
        end
    else
        for ψ in Ψ
            ψ(L) * ψ(U) < 0 && union!(roots, find_zero(ψ, (L, U)))
        end
    end
    sort!(roots)
    nodes = Vector{Float64}()
    weights = Vector{Float64}()
    for (l, r) in zip(roots[1:end-1], roots[2:end])
        val = [ψ((l + r) / 2.0) for ψ in Ψ]
        if all((val .* signs) .≥ 0)
            x, w = gausslegendre(q)
            @. x = (x + 1.0) / 2.0 * (r - l) + l
            @. w = w / 2.0 * (r - l)
            append!(nodes, x)
            append!(weights, w)
        end
    end
    return nodes, weights
end

function quadratureNodesWeights(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, rec::HyperRectangle{D}, q, surf, ∇Ψ, level=0) where {D}
    @assert !surf || (D > 1)
    # FIXME: the maximum depth should be a parameter (and why should it depend
    # on D)
    if level ≥ 20*D
        @warn "Maximum depth reached: resorting to low-order quadrature" maxlog=1
        if surf
            return [center(rec)], prod(i->high_corner(rec)[i] - low_corner(rec)[i], D-1)
        else
            return [center(rec)], prod(high_corner(rec).-low_corner(rec))
        end
    end
    ##### Pruning #####
    delInd = Vector{Int}()
    n = length(Ψ)
    for (i, ψ) in enumerate(Ψ)
        δ = bound(ψ, rec)
        if abs(ψ(center(rec))) ≥ δ
            if signs[i] * ψ(center(rec)) ≥ 0
                append!(delInd, i)
                n -= 1
            else
                return Vector{SVector{D,Float64}}(), Vector{Float64}()
            end
        end
    end
    deleteat!(signs, delInd)
    if n == 0
        if surf
            return Vector{SVector{D,Float64}}(), Vector{Float64}()
        else
            return tensorquad(rec, gausslegendre(q)...)
        end
    end
    deleteat!(Ψ, delInd)
    if ∇Ψ !== nothing
        deleteat!(∇Ψ, delInd)
    end
    ###################

    ##### Base case #####
    if D == 1
        return dim1NodesWeights(Ψ, signs, low_corner(rec)[1], high_corner(rec)[1], q)
    end
    ######################

    ##### Determine direction #####
    if ∇Ψ === nothing
        ∇Ψ = [∇(ψ, D) for ψ in Ψ]
    end
    k = argmax(abs.(∇Ψ[1](center(rec))))
    Ψ̃ = Vector{Function}()
    new_signs = Vector{Integer}()
    ∇Ψ̃ = Vector{Function}()
    for (ψ, s, ∇ψ) in zip(Ψ, signs, ∇Ψ)
        g = ∇ψ(center(rec))
        δ = [bound(x -> ∇ψ(x)[j], rec) for j = 1:D]
        if abs(g[k]) > δ[k] && sum((g .+ δ) .^ 2) / (g[k] - δ[k])^2 < 20
            ψL = x -> ψ(insert(x, k, low_corner(rec)[k]))
            sL = sgn(g[k], s, surf, -1)
            ∇ψL = x -> deleteat(∇ψ(insert(x, k, low_corner(rec)[k])), k)
            ψU = x -> ψ(insert(x, k, high_corner(rec)[k]))
            sU = sgn(g[k], s, surf, 1)
            ∇ψU = x -> deleteat(∇ψ(insert(x, k, high_corner(rec)[k])), k)
            append!(Ψ̃, [ψL, ψU])
            append!(new_signs, [sL, sU])
            append!(∇Ψ̃, [∇ψL, ∇ψU])
        else
            rec1, rec2 = split(rec)
            X1, W1 = quadratureNodesWeights(copy(Ψ), copy(signs), rec1, q, surf, copy(∇Ψ), level+1)
            X2, W2 = quadratureNodesWeights(copy(Ψ), copy(signs), rec2, q, surf, copy(∇Ψ), level+1)
            return (append!(X1, X2), append!(W1, W2))
        end
    end
    ###############################

    ##### Recursion in one less dimension #####
    X, W = quadratureNodesWeights(Ψ̃, new_signs, section(rec, k), q, false, ∇Ψ̃, level)
    nodes = Vector{SVector{D,Float64}}()
    weights = Vector{Float64}()
    for (x, w) in zip(X, W)
        if surf
            y = find_zero(y -> Ψ[1](insert(x, k, y)), (low_corner(rec)[k], high_corner(rec)[k]))
            x̃ = insert(x, k, y)
            ∇ϕ = ∇Ψ[1](x̃)
            push!(nodes, x̃)
            push!(weights, w * LinearAlgebra.norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(insert(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1NodesWeights(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], q, false)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end

function quadratureNodesWeights(Ψ::Vector{BernsteinPolynomial{D}}, signs::Vector{<:Integer}, q, surf, ∇Ψ, level=0) where{D}
    @assert !surf || (D > 1)
    rec = Ψ[1].domain
    if level ≥ D*20
        @warn "Maximum depth reached: resorting to low-order quadrature" maxlog=1
        if surf
            return [center(rec)], prod(i->high_corner(rec)[i] - low_corner(rec)[i], D-1)
        else
            return [center(rec)], prod(high_corner(rec).-low_corner(rec))
        end
    end
    ##### Pruning #####
    delInd = Vector{Int}()
    n = length(Ψ)
    for (i, ψ) in enumerate(Ψ)
        m, M = bound(ψ)
        if m ≥ 0
            if signs[i] ≥ 0
                append!(delInd, i)
                n -= 1
            else
                return Vector{SVector{D,Float64}}(), Vector{Float64}()
            end
        elseif M ≤ 0
            if signs[i] ≤ 0
                append!(delInd, i)
                n -= 1
            else
                return Vector{SVector{D,Float64}}(), Vector{Float64}()
            end
        end
    end
    deleteat!(signs, delInd)
    if n == 0
        if surf
            return Vector{SVector{D,Float64}}(), Vector{Float64}()
        else
            return tensorquad(rec, gausslegendre(q)...)
        end
    end
    deleteat!(Ψ, delInd)
    ###################

    ##### Base case #####
    if D == 1
        nodes, weights = dim1NodesWeights([x->ψ(x) for ψ in Ψ], signs, low_corner(rec)[1], high_corner(rec)[1], q)
        return nodes, weights
    end
    ######################

    ##### Determine direction #####
    K = Vector{Integer}()
    for j in 1:D
        non_zero = true
        for ∇ψ in ∇Ψ
            m, M = bound(∇ψ[j])
            if m * M ≤ 0
                non_zero = false
                break
            end
        end
        non_zero && push!(K, j)
    end
    if length(K) == 0
        split_ax = argmax(high_corner(rec) - low_corner(rec))
        Ψ1 = Vector{BernsteinPolynomial{D}}()
        Ψ2 = Vector{BernsteinPolynomial{D}}()
        for ψ in Ψ
            ψ1, ψ2 = split(ψ, split_ax)
            push!(Ψ1, ψ1); push!(Ψ2, ψ2)
        end
        ∇Ψ1 = [∇(ψ1) for ψ1 in Ψ1]; ∇Ψ2 = [∇(ψ2) for ψ2 in Ψ2]
        X1, W1 = quadratureNodesWeights(Ψ1, copy(signs), q, surf, ∇Ψ1, level+1)
        X2, W2 = quadratureNodesWeights(Ψ2, copy(signs), q, surf, ∇Ψ2, level+1)
        return (append!(X1, X2), append!(W1, W2))
    elseif length(K) == 1
        k = K[1]
    else
        c = svector(i->0.5, D); Grad = [[abs(∇ψ[j](c)) for j in 1:D] for ∇ψ in ∇Ψ]
        k = K[1]; s = sum(grad->grad[k]/sum(grad), Grad)
        for k̃ in K[2:end]
            s̃ = sum(grad->grad[k̃]/sum(grad), Grad)
            if s̃ > s
                k = k̃; s = s̃
            end
        end
    end
    Ψ̃ = Vector{BernsteinPolynomial{D-1}}()
    new_signs = Vector{Integer}()
    ∇Ψ̃ = Vector{SVector{D-1,BernsteinPolynomial{D-1}}}()
    for (ψ, s, ∇ψ) in zip(Ψ, signs, ∇Ψ)
        pos_neg = bound(∇ψ[k])[1] > 0 ? 1 : -1
        ψL = lower_restrict(ψ, k); sL = sgn(pos_neg, s, surf, -1); ∇ψL = svector(j->j<k ? lower_restrict(∇ψ[j], k) : lower_restrict(∇ψ[j+1], k), D-1)
        ψU = upper_restrict(ψ, k); sU = sgn(pos_neg, s, surf, 1);  ∇ψU = svector(j->j<k ? upper_restrict(∇ψ[j], k) : upper_restrict(∇ψ[j+1], k), D-1)
        append!(Ψ̃, [ψL, ψU])
        append!(new_signs, [sL, sU])
        append!(∇Ψ̃, [∇ψL, ∇ψU])
        # append!(∇Ψ̃, [∇(ψL), ∇(ψU)])
    end
    ###############################

    ##### Recursion in one less dimension #####
    X, W = quadratureNodesWeights(Ψ̃, new_signs, q, false, ∇Ψ̃, level)
    nodes = Vector{SVector{D,Float64}}()
    weights = Vector{Float64}()
    for (x, w) in zip(X, W)
        if surf
            y = find_zero(y -> Ψ[1](insert(x, k, y)), (low_corner(rec)[k], high_corner(rec)[k]))
            x̃ = insert(x, k, y)
            ∇ϕ = [∇Ψ[1][j](x̃) for j = 1:D]
            push!(nodes, x̃)
            push!(weights, w * LinearAlgebra.norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(insert(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1NodesWeights(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], q, false)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end
