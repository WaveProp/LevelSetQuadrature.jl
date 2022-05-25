using Test
using LevelSetQuadrature
using LinearAlgebra
using SpecialFunctions

using LevelSetQuadrature: svector

@testset "Hypersphere area/volume" begin
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
            U = HyperRectangle(1.1*svector(i->-r,n),1.1*svector(i->r,n))
            ϕ  = (x) -> sum(x .* x) - r^2
            ∇ϕ = (x) -> svector(i->2*x[i],n)
            # volume test
            X, W = quadgen(ϕ, ∇ϕ, U, -1;order=p)
            @test isapprox(sum(W),sphere_volume(r,n);atol)
            # area test
            X, W = quadgen(ϕ, ∇ϕ, U,0;order=p)
            @test isapprox(sum(W),sphere_area(r,n);atol)
        end
    end
end
