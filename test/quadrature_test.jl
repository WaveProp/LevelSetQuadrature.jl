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
U = HyperRectangle(SVector(0., 0.), SVector(1., 1.))
ϕ(x) = x[1]^2 + (x[2] + 0.5)^2 - 1
X, W = quadratureNodesWeights([ϕ], [-1], U, 10)

plot([x[1] for x in X], [x[2] for x in X], seriestype=:scatter)
θ = π/6 : π/300 : π/2
plot!(cos.(θ), sin.(θ).-0.5)

@test sum(W) ≈ (π / 3. - √3. / 4.)/2.
# @test sum(W) ≈ π / 3.
#############################################################


##### 3D #####
U = HyperRectangle(SVector(0., 0., 0.), SVector(1., 2., 0.3))
ϕ(x) = x[1]^2 + x[2]^2 + (x[3] + 0.5)^2 - 1
X, W = quadratureNodesWeights([ϕ], [-1], U, 5, true)

plot([x[1] for x in X], [x[2] for x in X], [x[3] for x in X], seriestype=:scatter)
# # θ = LinRange(0, π/3, 50)
# # t = LinRange(0, 2π, 100)
# # @. plot!(sin(θ)'*cos(t), sin(θ)*sin(t), cos(θ) - 0.5)

# @test sum(W) ≈ (π / 3. - π / 8.)/4.
# @test sum(W) ≈ π / 4.