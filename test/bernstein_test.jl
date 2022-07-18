using Test
using LevelSetQuadrature
using StaticArrays
using SpecialFunctions
using DynamicPolynomials

using LevelSetQuadrature: svector, low_corner, high_corner, center, gradient

U = HyperRectangle((0.,-1.,-1.), (2.,2.,3.))
L = low_corner(U); H = high_corner(U)
a = zeros(3,3,3)
a[1,1,1] = -1.; a[3,1,1] = 1.; a[1,3,1] = 2.; a[1,1,3] = 3.
ϕ = power2bernstein(a, U) # x1^2 + 2x2^2 + 3x3^2 - 1
N = 5 # number of tests of each category

@testset "Algebra tests" begin
    @testset "Evaluation tests" begin
        f = (x) -> x[1]^2 + 2x[2]^2 + 3x[3]^2 - 1
        for _ in 1:N
            x = @SVector rand(3)
            x = x .* (H - L) .+ L
            @test ϕ(x) ≈ f(x)
        end
    end

    @testset "Derivation tests" begin
        ∇ϕ = gradient(ϕ)
        ∇f = SVector(x->2x[1], x->4x[2], x->6x[3])
        for _ in 1:N
            x = @SVector rand(3)
            x = x .* (H - L) .+ L
            for i in 1:3
                @test ∇ϕ[i](x) ≈ ∇f[i](x)
            end
        end
    end

    @testset "Splitting tests" begin
        c = center(U)
        for d in 1:3
            ϕ1, ϕ2 = split(ϕ, d)
            for _ in 1:N
                x = svector(i->i == d ? c[i] : rand()*(H[i]-L[i])+L[i], 3)
                x̂ = deleteat(x, d)
                @test LevelSetQuadrature.upper_restrict(ϕ1, d)(x̂) ≈ LevelSetQuadrature.lower_restrict(ϕ2, d)(x̂) ≈ ϕ(x)
            end
        end
    end
end

@testset "Hypersphere area/volume" begin
    sphere_volume(r,n) = π^(n/2)/gamma(n/2+1)*r^n
    sphere_area(r,n)   = 2*π^(n/2)/gamma(n/2)*r^(n-1)
    r = 1.
    p = 5
    atol = 1e-3
    for n in 2:4
        @testset "dimension $n" begin
            U = HyperRectangle(1.0*svector(i->-r,n),1.0*svector(i->r,n))
            @testset "Polynomial interface" begin
                c = zeros(ntuple(i->3, n))
                c[1] = -1
                for i in 0:n-1
                    c[2*3^i+1] = 1
                end
                ϕ = LevelSetQuadrature.Polynomial(c)
                X, W = quadgen(ϕ,U,:negative; order=p)
                @test isapprox(sum(W),sphere_volume(r,n);atol)
                # area test
                X, W = quadgen(ϕ,U,:zero; order=p)
                @test isapprox(sum(W),sphere_area(r,n);atol)
            end
            @testset "Bernstein interface" begin
                c = zeros(ntuple(i->3, n))
                c[1] = -1
                for i in 0:n-1
                    c[2*3^i+1] = 1
                end
                ϕ = power2bernstein(c, U)
                X, W = quadgen(ϕ,:negative; order=p)
                @test isapprox(sum(W),sphere_volume(r,n);atol)
                # area test
                X, W = quadgen(ϕ,:zero; order=p)
                @test isapprox(sum(W),sphere_area(r,n);atol)
            end
            @testset "DynamicPolynomials interface" begin
                @polyvar x[1:n]
                ϕ = sum(i->i^2,x) - r^2
                # volume test
                X, W = quadgen(ϕ,U,:negative; order=p)
                @test isapprox(sum(W),sphere_volume(r,n);atol)
                # area test
                X, W = quadgen(ϕ,U,:zero; order=p)
                @test isapprox(sum(W),sphere_area(r,n);atol)
            end
        end
    end
end
