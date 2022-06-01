using SafeTestsets

#TODO: add tests for linearization functionality
# @safetestset "Linearization tests" begin include("linearization_test.jl") end
@safetestset "Function tests" begin include("quadrature_test.jl") end
@safetestset "Bernstein tests" begin include("bernstein_test.jl") end
