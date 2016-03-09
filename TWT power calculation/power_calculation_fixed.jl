module PowerCalculation

using DataFrames.readtable, ArrayViews.view

export calpower

const epsilon0 = 8.854e-12
const c = 299792458.0
const a = 2.032e-3
const b = a/2
const dx = 0.01e-3
const dz = dx
const dA = dx*dz
const nx = round(Int, 1.15*2/dx*1e-3) + 1
const nz = round(Int, 0.57*2/dz*1e-3) + 1
const ndata = 21

function sumpower(data::Matrix{Float64})
    output = 0.0
    for x in data
        output += c*epsilon0*dA*x^2
    end
    output
end

function checker!(data::Matrix{Real}, filtered::Matrix{Float64})
    j = 1
    for i in 1:size(data, 1)
        if abs(data[i, 1]) <= 1.15 && abs(data[i, 2]) <= 0.57
            filtered[j, :] = view(data, i, 4:6)
            j += 1
        end
    end
end

function calpower()
    powers = Array{Float64}(ndata)
    filteredtable = Array{Float64}(nx*nz, 3)
    table1 = readtable("J:/TWT/power2/1.txt", skipstart=3, header=false, separator=' ')
    arrtable = convert(Array, table1)
    checker!(arrtable, filteredtable)
    powers[1] = sumpower(filteredtable)
    for i in 2:ndata
        table = readtable("J:/TWT/power2/" * string(i) * ".txt", skipstart=3, header=false, separator=' ')
        arrtable = convert(Array, table)
        checker!(arrtable, filteredtable)
        powers[i] = sumpower(filteredtable)
    end
    powerAvg = mean(powers)
end

end