using LevelSetQuadrature
using Test
using StaticArrays

α = 1.0
β = SVector(1.,2.,3.)
δ = SVector(1.,2.,3.)
ϵ = 0.1

u = LevelSetQuadrature.Linearization(α,β,δ,ϵ)
v = LevelSetQuadrature.Linearization(α,β,δ,ϵ)

@test (u + v).ϵ == ϵ
