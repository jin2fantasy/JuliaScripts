## Calculate cutoff frequencies of modes in the circular waveguide
## Date:2015-06-22
## Author:wonjin
## Reference: David. M. Pozar. Microwave Enginnering. 4th edition. p.124

Radius = 5.5e-3 # diameter of circular waveguide, in [m]
const c = 299792458 # speed of light, in [m/s]

using Calculus
using Roots
# find the nth root of derivatives of bessel functions
function BesselJPrimeRoots(m::Int, n::Int)
    besseljm(x) = besselj(m, x)
    dbesseljm = derivative(besseljm)
    dbjzeros = fzeros(dbesseljm, 0.0, 50) # limits can be changed
                                       # if n is high
    if dbjzeros[1] == 0.0
        deleteat!(dbjzeros, 1)
    end

    testonemore(x::Float64) = dbesseljm(x)*dbesseljm(nextfloat(x)) > 0 &&
                                dbesseljm(x)*dbesseljm(prevfloat(x)) > 0

    while testonemore(dbjzeros[1])
        deleteat!(dbjzeros, 1)
    end

    dbjzeros[n]
end

function BesselJRoots(m::Int, n::Int)
    besseljm(x) = besselj(m, x)
    bjzeros = fzeros(besseljm, 0.0, 50)
    bjzeros[1] == 0 && (deleteat!(bjzeros, 1))
    testonemore(x::Float64) = besseljm(x)*besseljm(nextfloat(x)) > 0 &&
                            besseljm(x)*besseljm(prevfloat(x)) > 0
    while testonemore(bjzeros[1])
        deleteat!(bjzeros, 1)
    end

    bjzeros[n]
end

Radius = [5.5e-3]
open("cutoff_results.txt", "w") do f
    for r in Radius
        println(f, "radius is $(r*1e3) mm")
        for (m, n) in [(0,1), (0,2), (0,3), (0,4), (8,1)]
            println(f, "    cutoff for TE$(m),$(n) mode is ",
                    c*BesselJPrimeRoots(m, n)/(2pi*r)/1e9, " GHz")
        end
        for (m, n) in [(8,1)]
            println(f, "    cutoff for TM$(m),$(n) mode is ",
                    c*BesselJRoots(m, n)/(2pi*r)/1e9, " GHz")
        end
    end
end
