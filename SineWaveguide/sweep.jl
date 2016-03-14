# addprocs(8)
# using PyPlot
# require("D:\\JuliaScripts\\SineWaveguide\\test.jl")
# N = 2001
# f = linspace(190e9, 270e9, N)
# @everywhere β0 = linspace(0, pi/p, 19)
# @everywhere i = 1
# @everywhere detMNPQ1 = [x->detMNPQ(x, β0[j]) for j in 1:19]
# alldata = Array(Float64, N, 0)
# while i < 20
#     val = pmap(detMNPQ1[i], f)
#     alldata = [alldata val]
#     i += 1
# end
# plot(f, alldata, ".")
# plot(f, alldata[:,2])


addprocs(8)
using Plots
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\test2.jl")
N=4001
@everywhere detMNPQ2(x) = detMNPQ(x, 180/180*pi/p)
f = linspace(0, 300e9, N)
# f = linspace(209.1e9, 209.6e9, N)
val2 = pmap(detMNPQ2, f)
plot!(f, val2)

show()
