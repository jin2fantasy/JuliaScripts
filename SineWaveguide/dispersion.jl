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

const kzmn = (repmat(m*pi, 1, NN) ./
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
    bmpam
    dmpcm
    Rp
    Rm
    βn
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
    R_index = !map(isapprox, abs(repmat(βn, 1, length(m))),
                abs(repmat(m'*pi/p, length(n), 1)))
    Rp = (im * R_index .* repmat(βn, 1, length(m)) .*
        (1 - exp(repmat(im*p*βn, 1, length(m))) .* repmat((-1).^m', length(n), 1)) ./
        ((repmat((βn.^2), 1, length(m)) - repmat((m'*pi/p).^2, length(n), 1))))
    Rp[~R_index] = p/2
    Rp[:,1] = Rp[:,1] + ~R_index[:,1]*p/2
    Rm = (-im*R_index.*repmat(βn, 1, length(m)) .*
        (1 - exp(repmat(-im*p*βn, 1, length(m))) .* repmat((-1).^m', length(n), 1)) ./
        (repmat((βn.^2), 1, length(m)) - repmat((m'*pi/p).^2, length(n), 1)))
    Rm[~R_index] = p/2
    Rm[:,1] = Rm[:,1] + ~R_index[:,1]*p/2

    param1 = Params(lmn, Fl, Fpl, Gl, Gpl, Fν, Fpν, Gν, Gpν, 0, 0, 0, 0, 0)
    bmpam = -Fpl(NN, g+h) ./ Gpl(NN, g+h)
    dmpcm = Fpl(NN, g+h) ./ Gpl(NN, g+h)
    bmpam[isnan(bmpam)] = -1
    dmpcm[isnan(dmpcm)] = 1
    calbmpam_dmpcm!(bmpam, dmpcm, param1)
    Fl1 = Fl(1, g)
    Fpl1 = Fpl(1, g)
    Gl1 = Gl(1, g)
    Gpl1 = Gpl(1, g)
    M = (eye(length(n)) * p^2 .* repmat(ν', length(n), 1) .* repmat(Fpν', length(n), 1)
        - repmat(Fν, 1, length(n)) * lmn[1, 1] * (Fpl1[1] + bmpam[1] * Gpl1[1]) /
        (Fl1[1] + bmpam[1] * Gl1[1]) .* repmat(Rp[:,1].', length(n), 1) .*
        repmat(Rm[:,1], 1, length(n)))
    N = (eye(length(n)) * p^2 .* repmat(ν', length(n), 1) .* repmat(Gpν', length(n), 1)
        - repmat(Gν, 1, length(n)) * lmn[1, 1] * (Fpl1[1] + bmpam[1] * Gpl1[1]) /
        (Fl1[1] + bmpam[1] * Gl1[1]) .* repmat(Rp[:,1].', length(n), 1) .*
        repmat(Rm[:,1], 1, length(n)))
    P = (eye(length(n)) * p^2 .* repmat(ν', length(n), 1) .* repmat(Fpν', length(n), 1)
        + repmat(Fν, 1, length(n)) * lmn[1, 1] * (Fpl1[1] - dmpcm[1] * Gpl1[1]) /
        (-Fl1[1] + dmpcm[1] * Gl1[1]) .* repmat(Rp[:,1].', length(n), 1) .*
        repmat(Rm[:,1], 1, length(n)) .*
        exp(im * (repmat(βn', length(n), 1) - repmat(βn, 1, length(n))) * s))
    Q = (-eye(length(n)) * p^2 .* repmat(ν', length(n), 1) .* repmat(Gpν', length(n), 1)
        - repmat(Gν, 1, length(n)) * lmn[1, 1] * (Fpl1[1] - dmpcm[1] * Gpl1[1]) /
        (-Fl1[1] + dmpcm[1] * Gl1[1]) .* repmat(Rp[:,1].', length(n), 1) .*
        repmat(Rm[:,1], 1, length(n)) .*
        exp(im * (repmat(βn', length(n), 1) - repmat(βn, 1, length(n))) * s))

    param2 = Params(lmn[:, 1], Fl1, Fpl1, Gl1, Gpl1, Fν, Fpν, Gν, Gpν,
                    bmpam, dmpcm, Rp, Rm, βn)
    result = caldetMNPQ(M, N, P, Q, param2)
    real(result)
end

function calbmpam_dmpcm!(bmpam, dmpcm, pp::Params)
    @inbounds for kk in NN:-1:2
        yn1 = g + h*(kk-1)/NN
        dn = p - p/pi*acos((NN-2kk+2)/NN)
        dn1 = p - p/pi*acos((NN-2kk+4)/NN)
        lmn = slice(pp.lmn, :, kk)
        Fl = pp.Fl(kk, yn1)
        Fpl = pp.Fpl(kk, yn1)
        Gl = pp.Gl(kk, yn1)
        Gpl = pp.Gpl(kk, yn1)

        bmpam_index = (bmpam .== -1)
        dmpcm_index = (dmpcm .== -1)

        Ymn1n1a = (  ~bmpam_index * dn1/dn .*
                                    (Fl + bmpam .* Gl) ./
                                    (lmn .* (Fpl + bmpam .* Gpl))   )
        Ymn1n1a[isnan(Ymn1n1a)] = 0
        Ymn1n1a -= bmpam_index ./ lmn * dn1/dn
        Ymn1n1b = (  ~dmpcm_index * dn1/dn .*
                                    (Fl + dmpcm .* Gl) ./
                                    (lmn .* (Fpl + dmpcm .* Gpl))   )
        Ymn1n1b[isnan(Ymn1n1b)] = 0
        Ymn1n1b -= dmpcm_index ./ lmn * dn1/dn

        lmn = slice(pp.lmn, :, kk-1)
        Fl = pp.Fl(kk-1, yn1)
        Fpl = pp.Fpl(kk-1, yn1)
        Gl = pp.Gl(kk-1, yn1)
        Gpl = pp.Gpl(kk-1, yn1)

        bmpam[:] = (    (Ymn1n1a .* lmn .* Fpl - Fl) ./
                        (Gl - Ymn1n1a .* lmn .* Gpl)     )
        dmpcm[:] = (    (Ymn1n1b .* lmn .* Fpl - Fl) ./
                        (Gl - Ymn1n1b .* lmn .* Gpl)     )

    end
    nothing
end

function caldetMNPQ(M, N, P, Q, pp::Params)
    lm1 = pp.lmn
    Fν = pp.Fν
    Fpν = pp.Fpν
    Gν = pp.Gν
    Gpν = pp.Gpν
    Fl1 = pp.Fl
    Fpl1 = pp.Fpl
    Gl1 = pp.Gl
    Gpl1 = pp.Gpl
    bmpam = pp.bmpam
    dmpcm = pp.dmpcm
    Rp = pp.Rp
    Rm = pp.Rm
    βn = pp.βn
    @inbounds for kk in 2:length(m)
        M -= (repmat(Fν, 1, length(n)) * 2lm1[kk] *
        (Fpl1[kk] + bmpam[kk] * Gpl1[kk]) / (Fl1[kk] + bmpam[kk] * Gl1[kk]) .*
        repmat(Rp[:,kk].', length(n), 1) .* repmat(Rm[:,kk], 1, length(n)))
        N -= (repmat(Gν, 1, length(n)) * 2lm1[kk] *
        (Fpl1[kk] + bmpam[kk] * Gpl1[kk]) / (Fl1[kk] + bmpam[kk] * Gl1[kk]) .*
        repmat(Rp[:,kk].', length(n), 1) .* repmat(Rm[:,kk], 1, length(n)))
        P += (repmat(Fν, 1, length(n)) * 2lm1[kk] *
        (Fpl1[kk] - dmpcm[kk] * Gpl1[kk]) / (-Fl1[kk] + dmpcm[kk] * Gl1[kk]) .*
        repmat(Rp[:,kk].', length(n), 1) .* repmat(Rm[:,kk], 1, length(n)) .*
        exp(im * (repmat(βn', length(n), 1) - repmat(βn, 1, length(n))) * s))
        Q -= (repmat(Gν, 1, length(n)) * 2lm1[kk] *
        (Fpl1[kk] - dmpcm[kk] * Gpl1[kk]) / (-Fl1[kk] + dmpcm[kk] * Gl1[kk]) .*
        repmat(Rp[:,kk].', length(n), 1) .* repmat(Rm[:,kk], 1, length(n)) .*
        exp(im * (repmat(βn', length(n), 1) - repmat(βn, 1, length(n))) * s))
    end
    det([M N; P Q])
end
