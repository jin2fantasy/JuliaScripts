# constants
const c = 299792458
const μ0 = 4e-7pi
const ϵ0 = 1/(μ0*c^2)
# geometry parameters
const a = 0.77e-3
const p = 0.46e-3
const h = 0.43e-3
const g = (0.15e-3)/2
const s = p/2
# space harmonics number in region I
const n = -3:3
# standing wave number n region II and III
const m = 0:4
# number of steps in region II and III
const NN = 1000

const kzmn = (repmat(m, 1, NN) ./
repmat(p - p/pi*acos((NN-2(1:NN)'+2)/NN), length(m), 1))
const kx = pi/a

immutable Params
    lmn
    Fl
    Fpl
    Gl
    Gpl
    Fν
    Fpν
    Gν
    Gpν
end

function detMNPQ(freq, β0)
    k = 2pi*freq/c
    βn = β0 + 2pi*n/p

    ν = sqrt(abs(k^2 - βn.^2 - kx^2))
    ν_index = (k^2 - βn.^2 - kx^2) .< 0
    lmn = sqrt(abs(k^2 - kzmn.^2 - kx^2))
    lmn_index = (k^2 - kzmn.^2 - kx^2) .< 0

    Fν = ν_index .* sinh(ν * g) + ~ν_index .* sin(ν * g)
    Fpν = ν_index .* cosh(ν * g) + ~ν_index .* cos(ν * g)
    Gν = ν_index .* cosh(ν * g) + ~ν_index .* cos(ν * g)
    Gpν = ν_index .* sinh(ν * g) - ~ν_index .* sin(ν * g)
    Fl(nn, x) = (lmn_index[:, nn] .* sinh(lmn[:, nn] .* x) +
    ~lmn_index[:, nn] .* sin(lmn[:, nn] .* x))
    Fpl(nn, x) = (lmn_index[:, nn] .* cosh(lmn[:, nn] .* x) +
    ~lmn_index[:, nn] .* cos(lmn[:, nn] .* x))
    Gl(nn, x) = (lmn_index[:, nn] .* cosh(lmn[:, nn] .* x) +
    ~lmn_index[:, nn] .* cos(lmn[:, nn] .* x))
    Gpl(nn, x) = (lmn_index[:, nn] .* sinh(lmn[:, nn] .* x) -
    ~lmn_index[:, nn] .* sin(lmn[:, nn] .* x))

    R_index = abs(repmat(βn, 1, length(m))) .!= abs(repmat(m'*pi/p, length(n), 1))
    Rp = (im * R_index .* repmat(βn, 1, length(m)) .*
    (1 - (cos(repmat(βn, 1, length(m))*p) +
    im*sin(repmat(βn, 1, length(m))*p)) .* repmat((-1).^m', length(n), 1)) ./
    ((repmat((βn.^2), 1, length(m)) - repmat((m'*pi/p).^2, length(n), 1))))
    Rp[~R_index] = p/2
    Rp[:,1] = Rp[:,1] + ~R_index[:,1]*p/2
    Rm = -im*R_index.*repmat(βn, 1, length(m)) .*
    (1 - (cos(repmat(βn, 1, length(m))*p) -
    im*sin(repmat(βn, 1, length(m))*p)) .* repmat((-1).^m', length(n), 1)) ./
    ((repmat((βn.^2), 1, length(m)) - repmat((m'*pi/p).^2, length(n), 1)))
    Rm[~R_index] = p/2
    Rm[:,1] = Rm[:,1] + ~R_index[:,1]*p/2

    param1 = Params(lmn, Fl, Fpl, Gl, Gpl, Fν, Fpν, Gν, Gpν)
    bmpam = -Fpl(NN, g+h) ./ Gpl(NN, g+h)
    dmpcm = Fpl(NN, g+h) ./ Gpl(NN, g+h)
    bmpam[isnan(bmpam)] = -1;
    dmpcm[isnan(dmpcm)] = 1;
    calbmpam_dmpcm!(bmpam, dmpcm, param1)
    bmpam, dmpcm
end

function calbmpam_dmpcm!(bmpam, dmpcm, pp::Params)
    yn1 = g + h*(NN-1)/NN
    dn = p - p/pi * acos((-NN+2)/NN)
    dn1 = p - p/pi * acos((-NN+4)/NN)
    lmn = slice(pp.lmn, :, NN)
    Fl = pp.Fl(NN, yn1)
    Fpl = pp.Fpl(NN, yn1)
    Gl = pp.Gl(NN, yn1)
    Gpl = pp.Gpl(NN, yn1)

    for kk in NN:-1:2
        bmpam_index = (bmpam .== -1)
        Ymn1n1 = (  ~bmpam_index * dn1/dn .*
                                    (Fl + bmpam .* Gl) ./
                                    (lmn .* (Fpl + bmpam .* Gpl))   )
        Ymn1n1[isnan(Ymn1n1)] = 0;
        Ymn1n1 -= bmpam_index ./ lmn * dn1/dn;

#=        lmn = slice(pp.lmn, :, kk-1)
        Fl = pp.Fl(kk-1, yn1)
        Fpl = pp.Fpl(kk-1, yn1)
        Gl = pp.Gl(kk-1, yn1)
        Gpl = pp.Gpl(kk-1, yn1) =#

        bmpam[:] = (    (Ymn1n1 .* lmn .* Fpl - Fl) ./
                        (Gl - Ymn1n1 .* lmn .* Gpl)     )
        Ymn1n1 = (      dn1/dn * ((Fl + dmpcm .* Gl) ./
                                (lmn .* (Fpl + dmpcm .* Gpl)))      )
        dmpcm[:] = (    (Ymn1n1 .* lmn .* Fpl - Fl) ./
                        (Gl - Ymn1n1 .* lmn .* Gpl)     )

        yn1 -= h/NN
        dn = dn1
        dn1 = p - p/pi * acos((NN-2*kk+4)/NN)
        Fl = pp.Fl(kk-1, yn1)
        Fpl = pp.Fpl(kk-1, yn1)
        Gl = pp.Gl(kk-1, yn1)
        Gpl = pp.Gpl(kk-1, yn1)
    end
    nothing
end
