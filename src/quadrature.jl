function V(x, k, y)
    # FIXME: why this hack?
    x = isa(x, SVector) ? x : SVector(x)
    insert(x, k, y)
end

function Sgn(m, s, S, σ)
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

function tensorGaussQuadrature(rec::HyperRectangle{D}, GX, GW) where {D}
    if D == 1
        nodes = (GX .+ 1.0) ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1]) .+ low_corner(rec)[1]
        weights = GW ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1])
        return nodes, weights
    else
        rec1 = section(rec, 1)
        nodes = Vector{SVector{D,Float64}}()
        weights = Vector{Float64}()
        X, W = tensorGaussQuadrature(rec1, GX, GW)
        for (x, w) in zip(X, W)
            Y = (GX .+ 1.0) ./ 2.0 * (high_corner(rec)[1] - low_corner(rec)[1]) .+ low_corner(rec)[1]
            Ω = GW ./ 2.0 .* (high_corner(rec)[1] - low_corner(rec)[1])
            for (y, ω) in zip(Y, Ω)
                push!(nodes, V(x, 1, y))
                push!(weights, w * ω)
            end
        end
    end
    nodes, weights
end

function dim1NodesWeights(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, L, U, q)
    roots = [L, U]
    for ψ in Ψ
        union!(roots, find_zeros(ψ, L, U))
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

function quadratureNodesWeights(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, rec::HyperRectangle{D}, q, surf=false, ∇Ψ=nothing) where {D}
    @assert !surf || (D > 1)
    # plotrec(rec)
    ##### Pruning #####
    delInd = Vector{Int}()
    n = length(Ψ) ## n = N doesn't work? ##
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
    if n == 0
        if !surf
            return tensorGaussQuadrature(rec, gausslegendre(q)...)
        else
            return Vector{SVector{D,Float64}}(), Vector{Float64}()
        end
    end
    deleteat!(Ψ, delInd)
    deleteat!(signs, delInd)
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
            ψL = x -> ψ(V(x, k, low_corner(rec)[k]))
            sL = Sgn(g[k], s, surf, -1)
            ∇ψL = x -> deleteat(∇ψ(V(x, k, low_corner(rec)[k])), k)
            ψU = x -> ψ(V(x, k, high_corner(rec)[k]))
            sU = Sgn(g[k], s, surf, 1)
            ∇ψU = x -> deleteat(∇ψ(V(x, k, high_corner(rec)[k])), k)
            append!(Ψ̃, [ψL, ψU])
            append!(new_signs, [sL, sU])
            append!(∇Ψ̃, [∇ψL, ∇ψU])
        else
            rec1, rec2 = split(rec)
            X1, W1 = quadratureNodesWeights(copy(Ψ), copy(signs), rec1, q, surf, copy(∇Ψ))
            X2, W2 = quadratureNodesWeights(copy(Ψ), copy(signs), rec2, q, surf, copy(∇Ψ))
            return (append!(X1, X2), append!(W1, W2))
        end
    end
    ###############################

    ##### Recursion in one less dimension #####
    X, W = quadratureNodesWeights(Ψ̃, new_signs, section(rec, k), q, false, ∇Ψ̃)
    nodes = Vector{SVector{D,Float64}}()
    weights = Vector{Float64}()
    for (x, w) in zip(X, W)
        if surf
            y = find_zero(y -> Ψ[1](V(x, k, y)), (low_corner(rec)[k], high_corner(rec)[k]))
            x̃ = V(x, k, y)
            ∇ϕ = ∇Ψ[1](x̃)
            push!(nodes, x̃)
            push!(weights, w * LinearAlgebra.norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(V(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1NodesWeights(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], q)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, V(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end
