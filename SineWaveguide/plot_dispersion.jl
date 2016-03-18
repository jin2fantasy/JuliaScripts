addprocs(8)
using Plots
using DataFrames
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\mysolver.jl")
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\test3.jl")
@everywhere N = 361
@everywhere β = linspace(0, 2pi/p, N)
@everywhere detMNPQ1 = [x -> detMNPQ(x, β[i]) for i in 1:N]
f = pmap(mysolver, detMNPQ1, repmat([280e9], N), repmat([380e9], N))
pyplot(size = (600,400))
plot(β*p/pi*180, 1e-9f, line=(:solid, 2, :black))
simdis = readtable("WR03-mode1.txt", header=false, skipstart=2, separator=' ')
plot!(0:10:360, simdis[2], line=:none, marker=:circle)
show()

for (index,item) in enumerate(f)
    if !(typeof(item) <: Number)
        println(index, "\n", item)
    end
end

open("dispersion.txt", "w") do file
    for x in f
        write(file, string(x))
        write(file, "\n")
    end
end
