"""
    tensorquad(rec::HyperRectangle,x,w)

Given the one-dimensional quadrature nodes `x` and  weights `w` for integrating
a function between `0` and `1`, return the `nodes` and `weights` of the
corresponding tensor-product quadrature on `rec`.
"""
function tensorquad(rec::HyperRectangle{D}, x1d, w1d) where {D}
    xl, xu = low_corner(rec), high_corner(rec)
    μ      = prod(xu-xl)
    nodes = map(Iterators.product(ntuple(i->x1d,D)...)) do x̂
        # map reference nodes to rec
        rec(SVector(x̂))
    end |> vec
    weights = map(Iterators.product(ntuple(i->w1d,D)...)) do ŵ
        # scale reference weights
        prod(ŵ)*μ
    end |> vec
    nodes, weights
end

# base case where we compute the one dimensional quadrature
function dim1quad(Ψ::Vector{<:Function}, signs::Vector{<:Integer}, L, U, x1d, w1d, multi_zeros=true)
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
    nodes   = Vector{Float64}()
    weights = Vector{Float64}()
    for (l, r) in zip(roots[1:end-1], roots[2:end])
        val = [ψ((l + r) / 2.0) for ψ in Ψ]
        if all((val .* signs) .≥ 0)
            # push then rescale the nodes
            append!(nodes, x1d)
            append!(weights, w1d)
            is = length(nodes)-length(x1d)+1 # where the new nodes start
            for i in is:length(nodes)
                nodes[i]   = nodes[i] * (r - l) + l
                weights[i] = weights[i] * (r - l)
            end
        end
    end
    return nodes, weights
end

struct Parameters
    maxdepth::Int
    maxslope::Float64
end

"""
    quadgen(ϕ,∇ϕ,U,s;order=5,maxdepth=20,maxgradient=20)

TODO: extensively document this function
"""
function quadgen(ϕ::Function,U::HyperRectangle{N},s,∇ϕ=gradient(ϕ,Val(N));kwargs...) where {N}
    if s == :interior
        signs = [-1]
    elseif s == :exterior
        signs = [1]
    elseif s == :surface
        signs = [0]
    else
        error("unrecognized argument $s. Options are `:exterior`, `:interior`, and `:surface`")
    end
    Ω = ImplicitDomain([ϕ],[∇ϕ],signs,U)
    quadgen(Ω;kwargs...)
end

function quadgen(ϕ::BernsteinPolynomial,s;kwargs...)
    Ω  = BernsteinDomain([ϕ],[s])
    quadgen(Ω;kwargs...)
end

function quadgen(Ω::AbstractDomain;order=5,maxdepth=20,maxslope=10)
    par     = Parameters(maxdepth,maxslope)
    x1d,w1d = gausslegendre(order)
    # normalize quadrature from [-1,1] to [0,1] interval
    x1d .= (x1d .+ 1) ./ 2
    w1d .= w1d ./ 2
    surf = all(iszero,Ω.signs)
    _quadgen(Ω,surf,x1d,w1d,0,par)
end

function _quadgen(Ω::AbstractDomain{N,T},surf,x1d, w1d, level, par::Parameters) where {N,T}
    rec = Ω.rec
    xc  = center(rec)
    Ψ   = Ω.Ψ
    ∇Ψ  = Ω.∇Ψ
    signs  = Ω.signs
    D   = ambient_dimension(rec)
    @assert !surf || (D > 1)
    if level ≥ par.maxdepth
        @warn "Maximum depth reached: resorting to low-order quadrature"
        if surf
            return [xc], [prod(i->high_corner(rec)[i] - low_corner(rec)[i], D-1)]
        else
            return [xc], [prod(high_corner(rec).-low_corner(rec))]
        end
    end

    ctype = cell_type(Ω)
    # check for limiting cases of empty or whole cells
    if ctype == empty_cell
        return Vector{SVector{D,T}}(), Vector{T}()
    elseif ctype == whole_cell
        if surf
            return SVector{D,T}[], T[]
        else
            return tensorquad(rec,x1d,w1d)
        end
    end

    # base case
    if D == 1
        return dim1quad(Ψ, signs, low_corner(rec)[1], high_corner(rec)[1], x1d, w1d, true)
    end

    # find a heigh direction such that all of ∇Ψ are (provably) bounded away
    # from zero.
    bnds    = map(∇Ψ) do ∇ψ
        ntuple(d->bound(∇ψ[d],rec),D)
    end
    isvalid = ntuple(D) do dim
        all(bnds) do bnd
            (prod(bnd[dim])≥0) &&
            (1/minimum(abs,bnd[dim]) < par.maxslope)
        end
    end
    if !any(isvalid) # no valid direction so split
        Ω1,Ω2 = split(Ω)
        X1, W1 = _quadgen(Ω1, surf, x1d, w1d, level+1, par)
        X2, W2 = _quadgen(Ω2, surf, x1d, w1d, level+1, par)
        return (append!(X1, X2), append!(W1, W2))
    end

    # If there is a valid direction, we go down on it. Chose the direction which
    # is the least steep overall by maximizing the minimum of the derivative on
    # direction k over all functions
    ∇Ψc = map(∇Ψ) do ∇ψ
        ntuple(D) do d
            ∇ψc = abs.(∇ψ[d](xc))
            ∇ψc/norm(∇ψc)
        end
    end
    k = argmax(1:D) do dim
        if isvalid[dim]
            minimum(∇ψc -> ∇ψc[dim],∇Ψc)
        else
            -Inf
        end
    end
    Ω̃  = restrict(Ω,k,surf) # the D-1 dimensional domain
    X, W = _quadgen(Ω̃, false, x1d, w1d, level, par)
    nodes = Vector{SVector{D,T}}()
    weights = Vector{T}()
    for (x, w) in zip(X, W)
        if surf
            # if we get here there should be only one Ψ
            @assert length(Ψ) == 1
            ψ = first(Ψ)
            ∇ψ = first(∇Ψ)
            y = find_zero(y -> ψ(insert(x, k, y)), (low_corner(rec)[k], high_corner(rec)[k]))
            x̃ = insert(x, k, y)
            ∇ϕ = map(f->f(x̃),∇ψ)
            push!(nodes, x̃)
            push!(weights, w * norm(∇ϕ) / abs(∇ϕ[k]))
        else
            Φ = [y -> ψ(insert(x, k, y)) for ψ in Ψ]
            Y, Ω = dim1quad(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], x1d, w1d, false)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end

function _quadgen(Ψ::Vector{BernsteinPolynomial{D,T}}, signs::Vector{<:Integer}, rec, surf, ∇Ψ, x1d, w1d, level, par::Parameters) where{D,T}
    @assert !surf || (D > 1)
    @assert rec == Ψ[1].domain
    if level ≥  par.maxdepth
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
            return tensorquad(rec, x1d, w1d)
        end
    end
    deleteat!(Ψ, delInd)
    ###################

    ##### Base case #####
    if D == 1
        nodes, weights = dim1quad([x->ψ(x) for ψ in Ψ], signs, low_corner(rec)[1], high_corner(rec)[1], x1d,w1d,true)
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
        rec1, rec2 = split(rec,split_ax)
        Ψ1 = empty(Ψ)
        Ψ2 = empty(Ψ)
        for ψ in Ψ
            ψ1, ψ2 = split(ψ, split_ax)
            push!(Ψ1, ψ1); push!(Ψ2, ψ2)
        end
        ∇Ψ1 = [gradient(ψ1) for ψ1 in Ψ1]; ∇Ψ2 = [gradient(ψ2) for ψ2 in Ψ2]
        X1, W1 = _quadgen(Ψ1, copy(signs), rec1, surf, ∇Ψ1, x1d, w1d, level+1, par)
        X2, W2 = _quadgen(Ψ2, copy(signs), rec2, surf, ∇Ψ2, x1d, w1d, level+1, par)
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
    Ψ̃ = Vector{BernsteinPolynomial{D-1,T}}()
    new_signs = Vector{Integer}()
    ∇Ψ̃ = Vector{SVector{D-1,BernsteinPolynomial{D-1,T}}}()
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
    X, W = _quadgen(Ψ̃, new_signs, section(rec, k), false, ∇Ψ̃, x1d, w1d, level, par)
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
            Y, Ω = dim1quad(Φ, signs, low_corner(rec)[k], high_corner(rec)[k], x1d, w1d, false)
            for (y, ω) in zip(Y, Ω)
                push!(nodes, insert(x, k, y))
                push!(weights, w * ω)
            end
        end
    end
    ###########################################
    return nodes, weights
end

# utility functions
function StaticArrays.insert(x̂::SVector{<:Any,<:LinearizationDual{N,T}},k,y::T) where {N,T}
    rec = first(x̂) |> domain
    ŷ = LinearizationDual(y,zero(SVector{N,T}),zero(T),rec)
    insert(x̂,k,ŷ)
end

function StaticArrays.insert(x::Real,k,y::Real)
    insert(SVector(x),k,y)
end

"""
    sgn(m,s,surface,side)

Compute the sign of the upper and lower restrictions of a function `ψ` with sign
`s`. The `side` variable is eithe `1` for the upper restriction, or `-1` for the
lower restriction, and `m` is the sign of `∇ψ ⋅ eₖ`.
"""
function sgn(m, s, surface, side)
    if m*s*side > 0 || surface
        if m * side > 0
            return 1
        else
            return -1
        end
    else
        return 0
    end
end
