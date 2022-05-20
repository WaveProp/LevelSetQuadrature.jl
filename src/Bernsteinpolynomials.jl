##### Bernstein Polynomials ##### 
struct BernsteinPolynomial{D}
    coeffs::Array{<:Real,D}
    degree::NTuple{D,Integer}
end

BernsteinPolynomial(c::Array{<:Real,D}) where{D} = BernsteinPolynomial(c, size(c).-1)

function BernsteinBase(k::NTuple{D,Integer}, i::CartesianIndex{D}) where{D}
    c = prod(j -> binomial(k[j], i[j]-1), 1:D)
    x = prod(j -> "(1-x_$j)^$(i[j]-1) (x_$j)^$(k[j]-i[j]+1) ", 1:D)
    string(c, x)
end

function BernsteinHTML(k::NTuple{D,Integer}, i::CartesianIndex{D}) where{D}
    c = prod(j -> binomial(k[j], i[j]-1), 1:D)
    x = prod(1:D) do j
        k[j] == 0 && return ""
        p = i[j] == 2 ? "" : string(i[j] - 1)
        q = i[j] == k[j] ? "" : string(k[j]-i[j]+1)
        i[j] == 1        && return "(1-<i>x</i><sub><i>$j</i></sub>)<sup><i>$q</i></sup>"
        i[j] == k[j] + 1 && return "<i>x</i><sub><i>$j</i></sub><sup><i>$p</i></sup>"
        "<i>x</i><sub><i>$j</i></sub><sup><i>$p</i></sup>(1-<i>x</i><sub><i>$j</i></sub>)<sup><i>$q</i></sup>"
    end
    (c, x)
end

function (p::BernsteinPolynomial{D})(d::Integer, x::Real) where{D}
    @assert 1 ≤ d ≤ D
    k = p.degree[d]
    k == 0 && return BernsteinPolynomial(selectdim(p.coeffs, d, 1) .* 1, ntuple(i->i<d ? i : i+1, D-1))
    n = size(p.coeffs)[d]
    B̃ = mapslices(p.coeffs; dims=d) do b
        x == 0 && return b[1]
        x == 1 && return k ≥ n ? 0 : b[n]
        for i in k:-1:1
            if i ≥ n
                @. b[1:n-1] = b[1:n-1]*(1-x) + b[2:n]*x
                b[n] = b[n]*(1-x)
            else
                @. b[1:i] = b[1:i]*(1-x) + b[2:i+1]*x      
            end
        end
        b[1]
    end
    D == 1 && return B̃[1]
    BernsteinPolynomial(selectdim(B̃, d, 1) .* 1, ntuple(i->i<d ? p.degree[i] : p.degree[i+1], D-1))
end

(p::BernsteinPolynomial{1})(x::Real) = p(1, x)

function (p::BernsteinPolynomial{D})(x::SVector{D}) where{D}
    for d in 1:D 
        p = p(1, x[d])
    end
    p
end

(P::SVector{N,BernsteinPolynomial})(d::Integer, x::Real) where{N} = svector(i->P[i](d, x), N)
(P::SVector{N,BernsteinPolynomial})(x::SVector) where{N} = svector(i->P[i](x), N)

function ∇(p::BernsteinPolynomial{D}) where{D}
    ntuple(D) do d
        n = size(p.coeffs)[d]
        k = p.degree[d]
        coeffs = mapslices(p.coeffs, dims=d) do b
            c = (b[2:n] .- b[1:n-1]) .* k
            if k ≥ n
                push!(c, -b[n]*k)                
            end
            c
        end
        BernsteinPolynomial(coeffs, ntuple(i->i==d ? p.degree[i]-1 : p.degree[i], D))
    end |> SVector
end

Base.show(io::IO, ::MIME"text/plain", p::BernsteinPolynomial) = println(io, join([BernsteinBase(p.degree, i) for i in CartesianIndices(p.coeffs)], "+ "))

function Base.show(io::IO, ::MIME"text/html", p::BernsteinPolynomial)
    str = Vector{String}()
    for i in CartesianIndices(p.coeffs)
        p.coeffs[i] == 0 && continue
        (c, x) = BernsteinHTML(p.degree, i)
        term = c * p.coeffs[i] == 1 && x ≠ "" ? x : string(c*p.coeffs[i], x)
        push!(str, term)
    end
    println(io, join(str, " + "))
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

function Base.split(p::BernsteinPolynomial{D}, d::Integer, α=0.5) where{D}
    @assert 1 ≤ d ≤ D
    k = p.degree[d]
    k == 0 && return p, p
    n = size(p.coeffs)[d]
    coeffs = mapslices(p.coeffs, dims=d) do b
        c1 = Vector{Real}(); c2 = Vector{Real}()
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
    p1 = BernsteinPolynomial(selectdim(coeffs, d, 1:k+1)*1, p.degree)
    p2 = BernsteinPolynomial(selectdim(coeffs, d, k+1:k+n)*1, p.degree)
    p1, p2
end

##### Conversion #####
function renormalize(a::Vector{<:Real}, l::Real, r::Real)
    n = length(a)
    ã = copy(a)
    for i in 0:n-2
        ã[n-i:n] .*= (r-l)
        ã[n-i-1:n-1] .+= ã[n-i:n] ./ (r-l) .* l 
    end
    ã
end

function renormalize(A::Array{<:Real,D}, rec::HyperRectangle{D}) where(D)
    L = low_corner(rec)
    R = high_corner(rec)
    Ã = copy(A)
    for d in 1:D 
        Ã = mapslices(Ã, dims=d) do a
            renormalize(a, L[d], R[d])
        end
    end
    Ã
end

function power2Berstein(a::Array{<:Real,D}, k=size(a)) where{D}
    b = zeros(k)
    for i in CartesianIndices(a)
        temp = zeros(Tuple([i[j] for j = 1:D]))
        for l in CartesianIndices(temp) 
            temp[l] = a[l] * prod(1:D) do j
                binomial(i[j]-1, l[j]-1) / binomial(k[j]-1, l[j]-1)
            end
        end
        b[i] = sum(temp)
    end
    BernsteinPolynomial(b)
end

export BernsteinPolynomial, ∇, renormalize, power2Berstein