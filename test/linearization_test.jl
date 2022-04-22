using LevelSetQuadrature
using Test
using StaticArrays

α = 1.0
β = SVector(1.,2.,3.)
δ = SVector(1.,2.,3.)
ϵ = 0.1

rec = HyperRectangle(SVector(0.,0.,0.), SVector(1.,1.,1.))

f(x) = x[1]*x[2]*x[3]

LevelSetQuadrature.bound(f, rec)