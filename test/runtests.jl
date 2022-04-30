using SafeTestsets

#TODO: add tests for linearization functionality
@safetestset "Linearization tests" begin include("linearization_test.jl") end

@safetestset "Quadrature tests" begin include("quadrature_test.jl") end
