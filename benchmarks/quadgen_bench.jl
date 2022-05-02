using LevelSetQuadrature
using BenchmarkTools

# 3d quadrature of sphere
p = 5
U = HyperRectangle(1.1*svector(i->-1,3),1.1*svector(i->1,3))
ϕ(x)  = sum(x .* x) - 1^2
∇ϕ(x) = svector(i->2*x[i],3)
# volume test

@btime quadratureNodesWeights($([ϕ]), [-1], $U, $p, false, $([∇ϕ]))
@btime quadratureNodesWeights($([ϕ]), [-1], $U, $p, true, $([∇ϕ]))
