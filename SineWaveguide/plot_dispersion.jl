addprocs(8)
using PyPlot
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\mysolver.jl")
@everywhere include("D:\\JuliaScripts\\SineWaveguide\\test.jl")
@everywhere N = 361
@everywhere β = linspace(0, 2pi/p, N)
@everywhere detMNPQ1 = [x -> detMNPQ(x, β[i]) for i in 1:N]
f = pmap(mysolver, detMNPQ1, repmat([190e9], N), repmat([280e9], N))
plot(β*p/pi*180, f, ".-")
grid()

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
