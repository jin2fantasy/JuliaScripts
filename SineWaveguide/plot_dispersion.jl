addprocs(8)
using PyPlot
@everywhere using Roots
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\test.jl")
@everywhere N = 201
@everywhere β = linspace(0, pi/p, N)
@everywhere f = similar(β)
@everywhere detMNPQ1 = [x -> detMNPQ(x, β[y]) for y in 1:N]
# initial guesses
ig = Float64[195e9 + (270e9 - 195e9)/N*i for i in 1:N]
bks = Vector{Float64}[[ig[i] - 25e9, ig[i] + 25e9] for i in 1:N]
@everywhere fzero_maxeval50(a, b, c) = fzero(a, b, c, maxeval=50)
f = pmap(fzero_maxeval50, detMNPQ1, ig, bks)
plot(β*p/pi*180, f, ".-")
grid()
