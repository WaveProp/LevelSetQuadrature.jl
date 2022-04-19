using Roots
using ForwardDiff
using StaticArrays
using LinearAlgebra
using FastGaussQuadrature

function V(x, k, y)
    x = isa(x, SVector) ? x : SVector(x)
    insert(x, k, y)
end

function Sgn(m, s, S, σ)
    if m*s*σ > 0 || S
        if m * σ > 0
            return 1
        else
            return -1
        end
    else
        return 0
    end
end

function quadratureNodesWeights(Ψ::Array{<:Function,N}, Sign::Array{<:Integer,N}, rec::HyperRectangle{D}, q, Surf=false) where {N,D}
    @assert !Surf || (D > 1 && N == 1)
    ##### Base case #####
    if D == 1
        roots = [rec.L[1], rec.U[1]]
        for ψ in Ψ
            append!(roots, find_zeros(ψ, rec.L[1], rec.U[1]))
        end
        sort!(roots)
        nodes   = Vector{Float64}()
        weights = Vector{Float64}()
        for (l, r) in zip(roots[1:end-1], roots[2:end])
            val = [ψ((l+r)/2.) for ψ in Ψ]
            if all((val .* Sign) .≥ 0)
                x, w = gausslegendre(q)
                @. x = (x + 1.)/2. * (r - l) + l
                @. w = w /2. * (r - l)
                append!(nodes, x)
                append!(weights, w)
            end
        end
        return nodes, weights
    end
    ######################

    ##### Pruning #####
    # delInd = Vector{Int}()
    # for (i, ψ) in enumerate(Ψ)
    #     δ = bound(ψ, rec)
    #     if abs(ψ(center(rec))) ≥ δ
    #         if Sign[i] * ψ(center(rec)) ≥ 0
    #             append!(delInd, i)
    #             global N -= 1
    #         else
    #             return Vector{SVector{D, Float64}}(), Vector{Float64}()
    #         end
    #     end
    # end
    # if N == 0
    #     ## tensor Gauss quadrature
    # end
    ###################

    ##### Determine direction #####
    k = argmax(abs.(ForwardDiff.gradient(Ψ[1], center(rec))))
    Ψ̃ = Vector{Function}()
    Sigñ = Vector{Integer}()
    for (ψ, s) in zip(Ψ, Sign)
        ∇ψ(x) = ForwardDiff.gradient(ψ, x)
        g = ∇ψ(center(rec))
        δ = [bound(x -> ∇ψ(x)[j], rec) for j = 1:D]
        if abs(g[k]) > δ[k] && sum((g .+ δ).^2)/(g[k] - δ[k])^2 < 20
            ψL = x -> ψ(V(x, k, rec.L[k])); sL = Sgn(g[k], s, Surf, -1)
            ψU = x -> ψ(V(x, k, rec.U[k])); sU = Sgn(g[k], s, Surf,  1)
            append!(Ψ̃, [ψL, ψU]); append!(Sigñ, [sL, sU])
        else
            rec1, rec2 = divide(rec)
            X1, W1 = quadratureNodesWeights(Ψ, Sign, rec1, q, Surf)
            X2, W2 = quadratureNodesWeights(Ψ, Sign, rec2, q, Surf)
            return(append!(X1, X2), append!(W1, W2))
        end
    end
    ###############################

    ##### Recursion in one less dimension #####
    X, W = quadratureNodesWeights(Ψ̃, Sigñ, section(rec, k), q)
    nodes   = Vector{SVector{D, Float64}}()
    weights = Vector{Float64}()
    for (x, w) in zip(X, W)
        if Surf
            y = find_zero(y -> Ψ[1](V(x, k, y)), (rec.L[k], rec.U[k]))
            x̃ = V(x, k, y)
            ∇ϕ = ForwardDiff.gradient(Ψ[1], x̃)
            push!(nodes, x̃)
            push!(weights, w * LinearAlgebra.norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(V(x, k, y)) for ψ in Ψ]
            Y, Ω = quadratureNodesWeights(Φ, Sign, HyperRectangle(SVector(rec.L[k]), SVector(rec.U[k])), q)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, V(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################

    return nodes, weights
end
