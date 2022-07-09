using Test
using LevelSetQuadrature
using LinearAlgebra
using SpecialFunctions

using LevelSetQuadrature: svector

@testset "Hypersphere" begin
    sphere_volume(r,n) = π^(n/2)/gamma(n/2+1)*r^n
    sphere_area(r,n)   = 2*π^(n/2)/gamma(n/2)*r^(n-1)
    r = 1.5
    p = 5
    atol = 1e-3
    for n in 2:4
        @testset "dimension $n" begin
            # FIXME: there is an error when the bounding box touches the
            # boundary of the surface, so it has to be taken slighly larger for
            # the moment
            U = HyperRectangle(1.2*svector(i->-r,n),1.1*svector(i->r,n))
            ϕ  = (x) -> sum(x .* x) - r^2
            # volume test
            @testset "volume" begin
                X, W = quadgen(ϕ, U,:negative;order=p)
                @test isapprox(sum(W),sphere_volume(r,n);atol)
            end
            @testset "surface" begin
                # area test
                X, W = quadgen(ϕ, U,:zero;order=p)
                @test isapprox(sum(W),sphere_area(r,n);atol)
            end
        end
    end
end
