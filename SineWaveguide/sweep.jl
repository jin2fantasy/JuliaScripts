addprocs(8)
using PyPlot
require("D:\\JuliaScripts\\SineWaveguide\\dispersion.jl")
N = 2001
f = linspace(190e9, 270e9, N)
@everywhere detMNPQ1(x) = detMNPQ(x, 30/180*pi/p)
val = pmap(detMNPQ1, f)
plot(f, val)
