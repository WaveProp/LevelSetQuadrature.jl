using Test
using LevelSetQuadrature
using LinearAlgebra
using StaticArrays
using SpecialFunctions

using LevelSetQuadrature: svector, V

@testset "Hypersphere area/volume" begin
    sphere_volume(r,n) = π^(n/2)/gamma(n/2+1)*r^n
    sphere_area(r,n)   = 2*π^(n/2)/gamma(n/2)*r^(n-1)
    r = 1.
    p = 5
    atol = 1e-3
    for n in 2:4
        @testset "dimension $n" begin
            U = HyperRectangle(1.0*svector(i->-r,n),1.0*svector(i->r,n))
            c = zeros(ntuple(i->3, n))
            c[1] = -1
            for i in 0:n-1
                c[2*3^i+1] = 1
            end
            ϕ = power2Berstein(c)
            # volume test
            X, W = quadratureNodesWeights([ϕ], [-1], U, p, false, [∇(ϕ)])
            @test isapprox(sum(W),sphere_volume(r,n);atol)
            # FIXME area test
            # X, W = quadratureNodesWeights([ϕ], [-1], U, p, true, [∇(ϕ)])
            # @test isapprox(sum(W),sphere_area(r,n);atol)
        end
    end
end