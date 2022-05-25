using LevelSetQuadrature
using BenchmarkTools
using LevelSetQuadrature: svector

# 3d quadrature of sphere
p = 10
U = HyperRectangle(1.1*svector(i->-1,3),1.1*svector(i->1,3))
ϕ  = (x) -> sum(x .* x) - 1^2
∇ϕ = (x) -> svector(i->2*x[i],3)
# volume quad
@btime quadgen($ϕ,$∇ϕ,$U,-1;order=$p) # 32 ms
# surface quad
@btime quadgen($ϕ,$∇ϕ,$U,0;order=$p) # 15ms

##
n = 3
c = zeros(ntuple(i->3, n))
c[1] = -1
for i in 0:n-1
    c[2*3^i+1] = 1
end
ϕ = power2Berstein(c, U)
# volume quad
@btime x,w = quadgen($ϕ,-1;order=$p)
# surface quad
@btime x,w = quadgen($ϕ,0;order=$p)
