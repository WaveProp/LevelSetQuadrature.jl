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

function tensorGaussQuadrature(rec::HyperRectangle{D}, GX, GW) where{D}
    if D == 1
        nodes   = (GX .+ 1.) ./ 2. .* (rec.U[1] - rec.L[1]) .+ rec.L[1]
        weights = GW ./ 2. .* (rec.U[1] - rec.L[1])
        return nodes, weights
    else
        rec1     = section(rec, 1)
        nodes   = Vector{SVector{D, Float64}}()
        weights = Vector{Float64}()
        X, W = tensorGaussQuadrature(rec1, GX, GW)
        for (x, w) in zip(X, W)
            Y = (GX .+ 1.) ./ 2. * (rec.U[1] - rec.L[1]) .+ rec.L[1]
            Ω = GW ./ 2. .* (rec.U[1] - rec.L[1])
            for (y, ω) in zip(Y, Ω)
                push!(nodes, V(x, 1, y))
                push!(weights, w * ω)
            end
        end
    end
    nodes, weights
end

function dim1NodesWeights(Ψ::Vector{<:Function}, Sign::Vector{<:Integer}, L, U, q)
    roots = [L, U]
        for ψ in Ψ
            union!(roots, find_zeros(ψ, L, U))
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

function plotrec(rec::HyperRectangle{D}) where {D}
    if D == 2
       x = [rec.L[1], rec.U[1], rec.U[1], rec.L[1], rec.L[1]]
       y = [rec.L[2], rec.L[2], rec.U[2], rec.U[2], rec.L[2]]
       plot!(x, y, linecolor="black", legend=false)
    end
end

function quadratureNodesWeights(Ψ::Array{<:Function,N}, Sign::Array{<:Integer,N}, rec::HyperRectangle{D}, q, Surf=false, ∇Ψ=nothing) where {N,D}
    @assert !Surf || (D > 1 && N == 1)
    # plotrec(rec)
    ##### Pruning #####
    delInd = Vector{Int}()
    n = length(Ψ) ## n = N doesn't work? ##
    for (i, ψ) in enumerate(Ψ)
        δ = bound(ψ, rec)
        if abs(ψ(center(rec))) ≥ δ
            if Sign[i] * ψ(center(rec)) ≥ 0
                append!(delInd, i)
                n -= 1
            else
                return Vector{SVector{D, Float64}}(), Vector{Float64}()
            end
        end
    end
    if n == 0
        if !Surf
            return tensorGaussQuadrature(rec, gausslegendre(q)...)
        else
            return Vector{SVector{D, Float64}}(), Vector{Float64}()
        end
    end
    deleteat!(Ψ, delInd)
    deleteat!(Sign, delInd)
    if ∇Ψ !== nothing
        deleteat!(∇Ψ, delInd)
    end
    ###################

    ##### Base case #####
    if D == 1
        return dim1NodesWeights(Ψ, Sign, rec.L[1], rec.U[1], q)
    end
    ######################

    ##### Determine direction #####
    if ∇Ψ === nothing
        ∇Ψ = [∇(ψ, D) for ψ in Ψ]
    end
    k = argmax(abs.(∇Ψ[1](center(rec))))
    Ψ̃ = Vector{Function}()
    Sigñ = Vector{Integer}()
    ∇Ψ̃ = Vector{Function}()
    for (ψ, s, ∇ψ) in zip(Ψ, Sign, ∇Ψ)
        g = ∇ψ(center(rec))
        δ = [bound(x -> ∇ψ(x)[j], rec) for j = 1:D]
        if abs(g[k]) > δ[k] && sum((g .+ δ).^2)/(g[k] - δ[k])^2 < 20
            ψL = x -> ψ(V(x, k, rec.L[k])); sL = Sgn(g[k], s, Surf, -1); ∇ψL = x -> deleteat(∇ψ(V(x, k, rec.L[k])), k)
            ψU = x -> ψ(V(x, k, rec.U[k])); sU = Sgn(g[k], s, Surf,  1); ∇ψU = x -> deleteat(∇ψ(V(x, k, rec.U[k])), k)
            append!(Ψ̃, [ψL, ψU]); append!(Sigñ, [sL, sU]); append!(∇Ψ̃, [∇ψL, ∇ψU])
        else
            rec1, rec2 = divide(rec)
            X1, W1 = quadratureNodesWeights(copy(Ψ), copy(Sign), rec1, q, Surf, copy(∇Ψ))
            X2, W2 = quadratureNodesWeights(copy(Ψ), copy(Sign), rec2, q, Surf, copy(∇Ψ))
            return(append!(X1, X2), append!(W1, W2))
        end
    end
    ###############################

    ##### Recursion in one less dimension #####
    X, W = quadratureNodesWeights(Ψ̃, Sigñ, section(rec, k), q, false, ∇Ψ̃ )
    nodes   = Vector{SVector{D, Float64}}()
    weights = Vector{Float64}()
    for (x, w) in zip(X, W)
        if Surf
            y = find_zero(y -> Ψ[1](V(x, k, y)), (rec.L[k], rec.U[k]))
            x̃ = V(x, k, y)
            ∇ϕ = ∇Ψ[1](x̃)
            push!(nodes, x̃)
            push!(weights, w * LinearAlgebra.norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(V(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1NodesWeights(Φ, Sign, rec.L[k], rec.U[k], q)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, V(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end

# ϕ(x) = x[1]^2 + x[2]^2 + (x[3])^2 - 1
# ψL = x -> ϕ(V(x, 2, 0.75))
# U2 = HyperRectangle(SVector(0., 0.), SVector(0.375, 0.75))
# bound(ψL, U)
# ψL(center(U))

# ψLL = x -> ψL(V(x, 1, 0.375))
# U1 = HyperRectangle(SVector(0.), SVector(0.375))
# bound(ψLL, U)

# U3 = HyperRectangle(SVector(0.,0.,0.), SVector(1.,1.,1.))
# f = y -> ϕ(V(SVector(0.,0.), 2, y))
# f(2.)
# bound(f, U3)
