addprocs(8)
using Plots
using DataFrames
@everywhere include("mysolver.jl")
@everywhere N = 361
@everywhere include("test3.jl")
@everywhere β = linspace(0, 2pi/p, N)
@everywhere detMNPQ1 = [x -> detMNPQ(x, β[i]) for i in 1:N]
f = pmap(mysolver, detMNPQ1, repmat([280e9], N), repmat([340e9], N))


## When a = 0.53 mm, b = 0.5 mm, h = 0.32 mm, p = 0.44 mm
pyplot(size = (800,600))
ph = 720
simdis = readtable("mode1-2.txt", header=false, skipstart=2, separator=' ')
plot(β*p/pi*180 + ph, 1e-9f, linewidth=4, label="Theory")
plot!((0:10:360) + ph, simdis[2], marker=(4), label="Simulation")
# electron beam dispersion
# reference: http://www.sciencedirect.com/science/article/pii/S1007570411003364
phaseinrad = linspace(0,2pi*4,36*3)
e = 1.6e-19
m = 9.11e-31
V = 11e3; # beam voltage
A = 4*e^2*V^2/m^2/c^4
beta = sqrt((-A+sqrt(A^2+4*A))/2)
I = 100e-3; # beam current
J = I/(2*g*a*0.9)
k2 = phaseinrad/p
v0 = beta*c
n0 = J/e/v0
wp = sqrt(4*pi*n0*e^2/m)
w = k2*v0+wp*(1-beta^2)^(3/4)
plot!(phaseinrad*180/pi, w/(2e9*pi),linecolor=:black,label="BeamV=$(V/1000)kV")

title!("Dispersion curve")
xaxis!("Phase shift")
yaxis!("Frequency (GHz)",
    guidefont=Plots.Font("Helvetica",14,:hcenter,:vcenter,0.0,RGB{U8}(0.0,0.0,0.0)))
xlims!(ph, ph+360)
ylims!(270, 400)
xticks!(ph:60:(ph+360),
    tickfont=Plots.Font("Helvetica",11,:hcenter,:vcenter,0.0,RGB{U8}(0.0,0.0,0.0)))



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
