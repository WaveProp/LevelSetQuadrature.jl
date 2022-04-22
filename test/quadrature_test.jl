using LevelSetQuadrature

using Roots
using ForwardDiff
using StaticArrays
using LinearAlgebra
using FastGaussQuadrature

using Test
using Plots

##### 1D #####
# U = HyperRectangle(SVector(-2.), SVector(2.))
# ϕ(x) = x^2 - 1
# ψ(x) = sin(x)

# quadratureNodesWeights([ϕ, ψ], [-1, 1], U, 3)
###############################################

##### 2D #####
U = HyperRectangle(SVector(-1., -1.), SVector(1., 1.))
ϕ(x)  = x[1]^2 + (x[2])^2 - 1
∇ϕ(x) = SVector(2x[1], 2(x[2]))
plot(xlim=(-1.6,1.6), ylim=(-1.6, 1.6), aspect_ratio=1)
# X, W = quadratureNodesWeights([ϕ], [-1], U, 5, true)
X, W = quadratureNodesWeights([ϕ], [-1], U, 5, true, [∇ϕ])

plot!([x[1] for x in X], [x[2] for x in X], seriestype=:scatter)
θ = 0 : π/300 : π/2
plot!(cos.(θ), sin.(θ))

# @test sum(W) ≈ π / 4.
@test sum(W) ≈ 2π
#############################################################


##### 3D #####
U = HyperRectangle(SVector(-1.5, -1.5, -1.5), SVector(1.5, 1.5, 1.5))
ϕ(x) = x[1]^2 + x[2]^2 + (x[3])^2 - 1
∇ϕ(x) = SVector(2x[1], 2x[2], 2x[3])
X, W = quadratureNodesWeights([ϕ], [-1], U, 5, true, [∇ϕ])


plot([x[1] for x in X], [x[2] for x in X], [x[3] for x in X], seriestype=:scatter)
# # θ = LinRange(0, π/3, 50)
# # t = LinRange(0, 2π, 100)
# # @. plot!(sin(θ)'*cos(t), sin(θ)*sin(t), cos(θ) - 0.5)

# @test sum(W) ≈ 4π / 3
@test sum(W) ≈ 4π


##### 4D #####
U = HyperRectangle(SVector(-1.5, -1.5, -1.5, -1.5), SVector(1.5, 1.5, 1.5, 1.5))
ϕ(x) = x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 - 1
∇ϕ(x) = SVector(2x[1], 2x[2], 2x[3], 2x[4])
X, W = quadratureNodesWeights([ϕ], [-1], U, 5, false, [∇ϕ])

@test sum(W) ≈ π^2 / 2
@test sum(W) ≈ 2π^2