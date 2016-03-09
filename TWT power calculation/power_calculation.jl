using DataFrames

# All lengths are in mm.
epsilon0 = 8.854e-15
c = 299792458e3
a = 2.032
b = a/2
dx = 0.01
dz = dx
dA = dx*dz
nx = round(Int, 1.15*2/dx) + 1
nz = round(Int, 0.57*2/dz) + 1
input = zeros(nx*nz, 6)
output = zeros(nx*nz, 6)
powers = zeros(21)
ndata = 21


for t in 1:ndata
    efields = readtable("J:/TWT/power2/"*string(t)*".txt", skipstart=3, header=false, separator=' ')
    efields = convert(Array, efields)
    m = size(efields, 1)
    j = 1
    k = 1

tic()
    for i in 1:m
        if abs(efields[i, 1]) <= 1.15 && abs(efields[i, 2]) <= 0.57
            output[k, :] = efields[i, 1:6]*1e-3
            k += 1
        end
    end
toc()

    powerIn = 0
    powerOut = 0
    for i in 1:size(output, 1)
        # powerIn = powerIn + c*epsilon0*dA*(input[i,4]^2+input[i,5]^2+input[i,6]^2)
        powerOut += c*epsilon0*dA*(output[i,4]^2+output[i,5]^2+output[i,6]^2)
    end
    powers[t] = powerOut
end


powerAvg = mean(powers)
