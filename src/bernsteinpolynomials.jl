"""
    BernsteinPolynomial{D,T}

TODO: docstring
"""
struct BernsteinPolynomial{D,T}
    coeffs::Array{T,D}
    degree::NTuple{D,Integer}
    domain::HyperRectangle{D,T}
end

order(p::BernsteinPolynomial) = p.degree

"""
    reference_cube(::Val{D})

Return the `[0,1]ᴰ` reference domain as a `HyperRectangle{D,Float64}`.
"""
reference_cube(::Val{D}) where {D} = HyperRectangle(svector(i->0., D), svector(i->1., D))
reference_cube(D) = HyperRectangle(svector(i->0., D), svector(i->1., D))

BernsteinPolynomial(c::Array{<:Real,D}) where {D} = BernsteinPolynomial(c, size(c).-1, reference_cube(Val(D)))

function lower_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}
    @assert 1 ≤ d ≤ D
    BernsteinPolynomial(selectdim(p.coeffs, d, 1).*1, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
end

function upper_restrict(p::BernsteinPolynomial{D}, d::Integer) where{D}
    @assert 1 ≤ d ≤ D
    if p.degree[d] ≥ size(p.coeffs)[d]
        BernsteinPolynomial(reshape([0.], ntuple(i->1,D-1)), ntuple(i->1,D-1), section(p.domain, d))
    else
        BernsteinPolynomial(selectdim(p.coeffs, d, size(p.coeffs)[d]).*1, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
    end
end

function partial_application(p::BernsteinPolynomial{D,T},d::Integer, x::Real) where{D,T}
    @assert 1 ≤ d ≤ D
    l = low_corner(p.domain)[d]; r = high_corner(p.domain)[d]
    x = (x - l) / (r - l)
    sz = size(p.coeffs)
    sz′ =  ntuple(i -> i < d ? sz[i] : sz[i+1], D-1)
    len =  sz[d]
    out = similar(p.coeffs,sz′)
    for I in CartesianIndices(out)
        idxs = ntuple(i->i < d ? (I[i]:I[i]) : i==d ? (1:len) : (I[i-1]:I[i-1]),D)
        c̃    = view(p.coeffs,idxs...)
        out[I] = evaluate(SVector(x),c̃,Val{1}(),1,len)
    end
    BernsteinPolynomial(out, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1), section(p.domain, d))
end

function (p::BernsteinPolynomial{D})(x::SVector{D}) where{D}
    l = low_corner(p.domain); r = high_corner(p.domain)
    x₀ = (x - l) ./ (r - l)
    evaluate(x₀, p.coeffs, Val{D}(), 1, length(p.coeffs))
end
(p::BernsteinPolynomial)(x) = p(SVector(x))

(P::SVector{N,<:BernsteinPolynomial})(x) where {N} = svector(i->P[i](SVector(x)), N)

function gradient(p::BernsteinPolynomial{D}) where{D}
    svector(D) do d
        n = size(p.coeffs)[d]
        k = p.degree[d]
        l = low_corner(p.domain)[d]; u = high_corner(p.domain)[d]
        coeffs = mapslices(p.coeffs, dims=d) do b
            c = (b[2:n] .- b[1:n-1]) .* k ./ (u-l)
            if k ≥ n
                push!(c, -b[n]*k/(u-l))
            end
            c
        end
        BernsteinPolynomial(coeffs, ntuple(i->i==d ? p.degree[i]-1 : p.degree[i], D), p.domain)
    end
end

function Base.show(io::IO, p::BernsteinPolynomial)
    lb = low_corner(p.domain)
    ub = high_corner(p.domain)
    print(io, "Bernestein order ", order(p), " polynomial on ",
          '[', lb[1], ',', ub[1], ']')
    for i = 2:length(lb)
        print(io, " × [", lb[i], ',', ub[i], ']')
    end
end

function bound(p::BernsteinPolynomial)
    M = maximum(p.coeffs)
    m = minimum(p.coeffs)
    if M < 0 && prod(size(p.coeffs)) < prod(p.degree.+1)
        M = 0
    end
    if m > 0 && prod(size(p.coeffs)) < prod(p.degree.+1)
         m = 0
    end
    m, M
end

function Base.split(p::BernsteinPolynomial{D,T}, d::Integer, α=0.5) where {D,T}
    @assert 1 ≤ d ≤ D
    k = p.degree[d]
    k == 0 && return p, p
    n = size(p.coeffs)[d]
    coeffs = mapslices(p.coeffs, dims=d) do b
        c1 = Vector{T}(); c2 = Vector{T}()
        for i in k:-1:1
            if i ≥ n
                push!(c1, b[1])
                @. b[1:n-1] = b[1:n-1]*(1-α) + b[2:n]*α
                b[n] = b[n]*(1-α)
            else
                push!(c1, b[1])
                pushfirst!(c2, b[i+1])
                @. b[1:i] = b[1:i]*(1-α) + b[2:i+1]*α
            end
        end
        push!(c1, b[1])
        append!(c1, c2)
        c1
    end
    split_point = low_corner(p.domain)[d] + (high_corner(p.domain)[d] - low_corner(p.domain)[d])*α
    rec1, rec2 = split(p.domain, d, split_point)
    p1 = BernsteinPolynomial(collect(selectdim(coeffs, d, 1:k+1)), p.degree, rec1)
    p2 = BernsteinPolynomial(collect(selectdim(coeffs, d, k+1:k+n)), p.degree, rec2)
    p1, p2
end

function rebase(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:n-2
        ã[n-i:n] .*= (r-l)
        ã[n-i-1:n-1] .+= ã[n-i:n] ./ (r-l) .* l
    end
    ã
end

function rebase(A::Array{<:Real,D}, rec::HyperRectangle{D}) where(D)
    L = low_corner(rec)
    R = high_corner(rec)
    Ã = copy(A)
    for d in 1:D
        Ã = mapslices(Ã, dims=d) do a
            rebase(a, L[d], R[d])
        end
    end
    Ã
end

"""
    power2bernstein

TODO: document and write a `jldoctest`
"""
function power2bernstein(a::Array{<:Real,D}, U::HyperRectangle{D}=□(D), k=size(a).-1) where{D}
    b = zeros(k.+1)
    a = rebase(a, U)
    for i in CartesianIndices(a)
        temp = zeros(Tuple([i[j] for j = 1:D]))
        for l in CartesianIndices(temp)
            temp[l] = a[l] * prod(1:D) do j
                binomial(i[j]-1, l[j]-1) / binomial(k[j], l[j]-1)
            end
        end
        b[i] = sum(temp)
    end
    BernsteinPolynomial(b, k, U)
end

# Adapted from FastChebInterp
@fastmath function evaluate(x::SVector{N}, c::AbstractArray, ::Val{dim}, i1, len) where {N,dim}
    n = size(c,dim)
    @inbounds xd = x[dim]
    # idea taken from here https://personal.math.ubc.ca/~cass/graphics/text/www/pdf/a6.pdf
    if dim == 1
        s = 1-xd
        @inbounds P = c[i1]
        C = (n-1)*xd
        for k in 1:n-1
            @inbounds P = P*s + C*c[i1+k]
            C = C*(n-k-1)/(k+1)*xd
        end
        return P
    else
        Δi = len ÷ n # column-major stride of current dimension

        # we recurse downward on dim for cache locality,
        # since earlier dimensions are contiguous
        dim′ = Val{dim-1}()

        s = 1-xd
        P = evaluate(x, c, dim′, i1, Δi)
        C = (n-1)*xd
        for k in 1:n-1
            P = P*s + C*evaluate(x,c,dim′,i1+k*Δi,Δi)
            C = C*(n-k-1)/(k+1)*xd
        end
        return P
    end
end
