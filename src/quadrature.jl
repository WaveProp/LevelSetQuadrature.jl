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

# helper function to map {:negative,:positive,:zero} to the necessary data
# required by quadgen
function _symbol_to_signs_and_surf(s::Symbol)
    if s == :negative
        signs = [-1]; surf = false
    elseif s == :positive
        signs = [1];  surf = false
    elseif s == :zero
        signs = [0];  surf = true
    else
        error("unrecognized argument $s. Options are `:positive`, `:negative`, and `:zero`")
    end
    return signs,surf
end

### Function interface
function quadgen!(X, W, ϕ::Function,U::HyperRectangle{N},s,∇ϕ=gradient(ϕ,Val(N));kwargs...) where {N}
    signs, surf = _symbol_to_signs_and_surf(s)
    Ω = ImplicitDomain([ϕ],[∇ϕ],signs,U)
    quadgen!(X, W, Ω,surf;kwargs...)
end

"""
    quadgen(ϕ,∇ϕ,U,s;order=5,maxdepth=20,maxgradient=20)

Return the quadrature nodes `X` and weights `W` for integrating over the domain `V∩U`
where `V = {ϕ < 0}` if `s = :negative`
      `V = {ϕ > 0}` if `s = :positive`
      `V = {ϕ = 0}` if `s = :surface`
"""
function quadgen(ϕ::Function,U::HyperRectangle{N},s,∇ϕ=gradient(ϕ,Val(N));kwargs...) where {N}
    X = Vector{SVector{N,Float64}}()
    W = Vector{Float64}()
    quadgen!(X, W, ϕ, U, s, ∇ϕ;kwargs...)
end
######################

### Polynomial interface
function quadgen!(X, W, p::Polynomial{N}, U::HyperRectangle{N}, s;kwargs...) where {N}
    coeffs = p.coeffs
    ϕ = power2bernstein(coeffs, U)
    quadgen!(X, W, ϕ, s;kwargs...)
end

function quadgen(p::Polynomial{N}, U::HyperRectangle{N}, s;kwargs...) where {N}
    X = Vector{SVector{N,Float64}}()
    W = Vector{Float64}()
    quadgen!(X, W, p, U, s;kwargs...)
end
########################

### Bernstein Polynomial interface
function quadgen!(X, W, ϕ::BernsteinPolynomial{N},s;kwargs...) where{N}
    signs,surf = _symbol_to_signs_and_surf(s)
    Ω  = BernsteinDomain([ϕ],signs)
    quadgen!(X, W, Ω, surf;kwargs...)
end

function quadgen(ϕ::BernsteinPolynomial{N},s;kwargs...) where{N}
    X = Vector{SVector{N,Float64}}()
    W = Vector{Float64}()
    quadgen!(X, W, ϕ, s;kwargs...)
end
##################################

### AbstractDomain interface
function quadgen!(X, W, Ω::AbstractDomain, surf;order=5,maxdepth=20,maxslope=10)
    par     = Parameters(maxdepth,maxslope)
    x1d,w1d = gausslegendre(order)
    # normalize quadrature from [-1,1] to [0,1] interval
    x1d .= (x1d .+ 1) ./ 2
    w1d .= w1d ./ 2
    X̃, W̃ = _quadgen(Ω,surf,x1d,w1d,0,par)
    append!(X, X̃)
    append!(W, W̃)
    X, W
end

function quadgen(Ω::AbstractDomain{N},surf;order=5,maxdepth=20,maxslope=10) where{N}
    X = Vector{SVector{N,Float64}}()
    W = Vector{Float64}()
    quadgen!(X, W, Ω, surf;order, maxdepth, maxslope)
end
############################


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
            (prod(bnd[dim])>0) &&
            (sum(bd->maximum(abs,bd), bnd)/minimum(abs,bnd[dim]) < par.maxslope)
        end
    end
    if !any(isvalid) # no valid direction so split
        Ω1,Ω2 = split(Ω)
        X1, W1 = _quadgen(Ω1, surf, x1d, w1d, level+1, par)
        X2, W2 = _quadgen(Ω2, surf, x1d, w1d, level+1, par)
        return (append!(X1, X2), append!(W1, W2))
    end

    # If there is a valid direction, we go down on it. Choose the direction which
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
            minimum(∇ψc -> abs(∇ψc[dim]),∇Ψc)
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
            lk, rk = low_corner(rec)[k], high_corner(rec)[k]
            ψₖ(y) = ψ(insert(x, k, y))
            if ψₖ(lk) * ψₖ(rk) < 0
                y = find_zero(ψₖ, (lk, rk))
                x̃ = insert(x, k, y)
                ∇ϕ = map(f->f(x̃),∇ψ)
                push!(nodes, x̃)
                push!(weights, w * norm(∇ϕ) / abs(∇ϕ[k]))
            end
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
