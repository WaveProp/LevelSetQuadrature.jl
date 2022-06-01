# TODO: add linearization tests
using LevelSetQuadrature
using BenchmarkTools
using LevelSetQuadrature: svector, ImplicitDomain, bound, linearization
using StaticArrays
import LevelSetQuadrature as LSQ
using ForwardDiff

# 3d quadrature of sphere
p = 10
U = HyperRectangle(1.1*svector(i->-1,3),1.1*svector(i->1,3))
x₀ = @SVector rand(3)
ϕ  = (x) -> sum(x .* x) - 1^2
ϕ  = (x) -> x[1]^2 + x[2]^2 + x[3]^2 - 1^2

LSQ.gradient(ϕ,x₀)

dp = LSQ.gradient(ϕ,Val(3))



# ϕ  = (x) -> x[1]*x[1] + x[2]*x[2] + x[3]*x[3] - 1^2
∇ϕ = svector(i->(x) -> 2*x[i],3)

∇ϕ[3](x₀)



bound(ϕ,U)

bound(∇ϕ[1],U)

∇ϕ = LSQ.gradient(ϕ,Val(3))
bound(∇ϕ[2],U)


∇ϕ = (x) -> ForwardDiff.gradient(ϕ,x)
bound(∇ϕ,U)

# using Zygote
# ∇ϕ = (x...) -> gradient(ϕ,x...)[1][1]
# bound(∇ϕ,U)
# L = linearization(∇ϕ, U)
# b = bound(∇ϕ,U)

# Yota
using ForwardDiff

# see https://discourse.julialang.org/t/how-to-do-partial-derivatives/19869/9
using ForwardDiff: derivative, Dual
import LevelSetQuadrature as LSQ

∇ϕ = LSQ.gradient(ϕ,Val(3))

bound(∇ϕ[3],U)

L = linearization(∇ϕ, U)
L = linearization(pf, U)

L = linearization(∇ϕ, U)

b = bound(∇ϕ,U)


f = (x) -> SVector(x[1],x[2],x[3])
