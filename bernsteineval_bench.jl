using BenchmarkTools
using Test

@fastmath function bernstein(x,c::Vector,k::Integer=length(c)-1,b=similar(c))
    copy!(b,c)
    n = length(b)
    x == 0 && return b[1]
    x == 1 && return k ≥ n ? 0. : b[n]
    for i in k:-1:1
        if i ≥ n
            for j in 1:n-1
                @inbounds b[j] = b[j]*(1-x) + b[j+1]*x
            end
            b[n] = b[n]*(1-x)
        else
            for j in 1:i
                @inbounds b[j] = b[j]*(1-x) + b[j+1]*x
            end
        end
    end
    b[1]
end

@fastmath function bernstein_horner(t,y,n=length(y)-1)
    s = 1-t
    P = y[1]
    C = n*t
    for k in 1:n
        @inbounds P = P*s + C*y[k+1]
        C = C*(n-k)/(k+1)*t
    end
    return P
end

@fastmath function bernstein_old(x, c::Vector, k::Integer = length(c)-1)
    b = copy(c)
    n = length(b)
    x == 0 && return b[1]
    x == 1 && return k ≥ n ? 0. : b[n]
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


n = 10

c = rand(n)
b = copy(c)

x = rand()

v1 = bernstein(x,c)
v2 = bernstein_old(x,c)
v3 = bernstein_horner(x,c)

@test v1 ≈ v2 ≈ v3

b1 = @btime bernstein($x,$c,$(n-1),$b)
b2 = @btime bernstein_old($x,$c)
b3 = @btime bernstein_horner($x,$c)
