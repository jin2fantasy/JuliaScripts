addprocs(8)
using Plots
using DataFrames
@everywhere include("mysolver.jl")
@everywhere include("test3.jl")
@everywhere N = 361
@everywhere β = linspace(0, 2pi/p, N)
@everywhere detMNPQ1 = [x -> detMNPQ(x, β[i]) for i in 1:N]
f = pmap(mysolver, detMNPQ1, repmat([280e9], N), repmat([380e9], N))
gadfly(size = (400,300))
simdis = readtable("WR03-mode1.txt", header=false, skipstart=2, separator=' ')
plot(β*p/pi*180, 1e-9f, linewidth=4, label="Theory")
scatter!(0:10:360, simdis[2], markersize=4, label="Simulation")
title!("Dispersion curve")
xaxis!("Phase shift")
yaxis!("Frequency (GHZ)",
    guidefont=Plots.Font("Helvetica",14,:hcenter,:vcetner,0.0,RGB{U8}(0.0,0.0,0.0)))
xlims!(0, 360)
xticks!(0:60:360,
    tickfont=Plots.Font("Helvetica",11,:hcenter,:vcetner,0.0,RGB{U8}(0.0,0.0,0.0)))
gui()

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
