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
# # space harmonics number in region I
# const n = -3:3
# # standing wave number n region II and III
# const m = 0:4
# number of steps in region II and III
const NN = 1000
# wave number in x
const kx = pi/a
#####################
yn(n) = g + h*n/NN
dn(n) = p - p/pi*acos((NN-2n+2)/NN)
k(f) = 2pi*f/c
kzmn(m, n) = m*pi / dn(n)
betan(b0, n) = b0 + 2pi*n/p
nu(f, b0, n) = sqrt(abs(k(f)^2 - betan(b0, n)^2 - kx^2))
nu_index(f, b0, n) = (k(f)^2 - betan(b0, n)^2 - kx^2) < 0
lmn(f, m, n) = sqrt(abs(k(f)^2 - kzmn(m ,n)^2 - kx^2))
lmn_index(f, m, n) = (k(f)^2 - kzmn(m, n)^2 - kx^2) < 0
Fnu(f, b0, n) = nu_index(f, b0, n) * sinh(nu(f, b0, n) * g) +
    ~nu_index(f, b0, n) * sin(nu(f, b0, n) * g)
Fpnu(f, b0, n) = nu_index(f, b0, n) * cosh(nu(f, b0, n) * g) +
    ~nu_index(f, b0, n) * cos(nu(f, b0, n) * g)
Gnu(f, b0, n) = nu_index(f, b0, n) * cosh(nu(f, b0, n) * g) +
    ~nu_index(f, b0, n)* cos(nu(f, b0, n) * g)
Gpnu(f, b0, n) = nu_index(f, b0, n) * sinh(nu(f, b0, n) * g) -
    ~nu_index(f, b0, n) * sin(nu(f, b0, n) * g)
Fl(f, m, n, x) = (lmn_index(f, m, n) * sinh(lmn(f, m, n) * x) +
    ~lmn_index(f, m, n) * sin(lmn(f, m, n) * x))
Fpl(f, m, n, x) = (lmn_index(f, m, n) * cosh(lmn(f, m, n) * x) +
    ~lmn_index(f, m, n) * cos(lmn(f, m, n) * x))
Gl(f, m, n, x) = (lmn_index(f, m, n) * cosh(lmn(f, m, n) * x) +
    ~lmn_index(f, m, n) * cos(lmn(f, m, n) * x))
Gpl(f, m, n, x) = (lmn_index(f, m, n) * sinh(lmn(f, m, n) * x) -
    ~lmn_index(f, m, n) * sin(lmn(f, m, n) * x))

function Rp(b0, n, m)
    if isapprox(abs(betan(b0, n)), abs(m*pi/p))
        m == 0 && return p
        return convert(Complex, p/2)
    else
        return (im * betan(b0, n) * (1 - exp(im * betan(b0, n) * p) * (-1)^m) /
                (betan(b0, n)^2 - (m*pi/p)^2))
    end
end
function Rm(b0, n, m)
    if isapprox(abs(betan(b0, n)), abs(m*pi/p))
        m == 0 && return p
        return p/2
    else
        return (-im * betan(b0, n) * (1 - exp(-im * betan(b0, n) * p) * (-1)^m) /
                (betan(b0, n)^2 - (m*pi/p)^2))
    end
end

function evalbmpam_dmpcm(f, m)
    bmpam = -Fpl(f, m, NN, g+h) / Gpl(f, m, NN, g+h)
    dmpcm = Fpl(f, m, NN, g+h) / Gpl(f, m, NN, g+h)
    @inbounds for kk in NN:-1:2
        bmpam_index = (bmpam == -1)
        dmpcm_index = (dmpcm == -1)
        yn1 = yn(kk-1)
        dn1 = dn(kk-1)
        dn0 = dn(kk)
        lmk = lmn(f, m, kk)
        lmk1 = lmn(f, m, kk-1)
        _Fl(x) = Fl(f, m, x, yn1)
        _Fpl(x) = Fpl(f, m, x, yn1)
        _Gl(x) = Gl(f, m, x, yn1)
        _Gpl(x) = Gpl(f, m, x, yn1)
        if bmpam_index
            Ymn1n1a = - dn1 / dn0 / lmk
        else
            Ymn1n1a = (dn1/dn0 *
            (_Fl(kk) + bmpam * _Gl(kk)) /
            (lmk * (_Fpl(kk) + bmpam * _Gpl(kk)))   )
        end
        if dmpcm_index
            Ymn1n1b = - dn1 / dn0 / lmk
        else
            Ymn1n1b = (dn1/dn0 *
            (_Fl(kk) + dmpcm * _Gl(kk)) /
            (lmk * (_Fpl(kk) + dmpcm * _Gpl(kk)))   )
        end

        bmpam = ( (Ymn1n1a .* lmk1 * _Fpl(kk-1) - _Fl(kk-1)) ./
                (_Gl(kk-1) - Ymn1n1a * lmk1 * _Gpl(kk-1))     )
        dmpcm = ( (Ymn1n1b .* lmk1 * _Fpl(kk-1) - _Fl(kk-1)) ./
                (_Gl(kk-1) - Ymn1n1b * lmk1 * _Gpl(kk-1))     )
    end
    (bmpam, dmpcm)::Tuple{Float64, Float64}
end

function detMNPQ(freq, beta0)
    n = -3:3
    np = -3:3 # for clarity
    m = 0:4
    M = Array{Complex{Float64}}(length(n), length(n))
    N, P, Q = similar(M), similar(M), similar(M)

    @inbounds for ii in 1:length(np), jj in 1:length(n), kk in 1:length(m)
        bmpam, dmpcm = evalbmpam_dmpcm(freq, m[kk])
        _lmn = lmn(freq, m[kk], 1)
        _nu = nu(freq, beta0, np[ii])
        _Fnu = Fnu(freq, beta0, n[jj])
        _Fpnu = Fpnu(freq, beta0, np[ii])
        _Gnu = Gnu(freq, beta0, n[jj])
        _Gpnu = Gpnu(freq, beta0, np[ii])
        _Fl = Fl(freq, m[kk], 1, g)
        _Fpl = Fpl(freq, m[kk], 1, g)
        _Gl = Gl(freq, m[kk], 1, g)
        _Gpl = Gpl(freq, m[kk], 1, g)
        _Rp = Rp(beta0, np[ii], m[kk])
        _Rm = Rm(beta0, n[jj], m[kk])
        expbeta = exp(im * (betan(beta0, np[ii]) - betan(beta0, n[jj]) * s))

        M[jj,ii] = ( (ii == jj) * p^2 * _nu * _Fpnu -
        _Fnu * 2 / (1 + (m[kk] == 0)) * _lmn * (_Fpl +
        bmpam * _Gpl) / (_Fl + bmpam * _Gl) * _Rp * _Rm )

        N[jj,ii] = ( (ii == jj) * p^2 * _nu * _Gpnu -
        _Gnu * 2 / (1 + (m[kk] == 0)) * _lmn * (_Fpl +
        bmpam * _Gpl) / (_Fl + bmpam * _Gl) * _Rp * _Rm )

        P[jj,ii] = ( (ii == jj) * p^2 * _nu * _Fpnu +
        _Fnu * 2 / (1 + (m[kk] == 0)) * _lmn * (_Fpl - dmpcm * _Gpl) /
        (-_Fl + dmpcm * _Gl) * _Rp * _Rm * expbeta )

        Q[jj,ii] = ( -(ii == jj) * p^2 * _nu * _Gpnu -
        _Gnu * 2 / (1 + (m[kk] == 0)) * _lmn * (_Fpl - dmpcm * _Gpl) /
        (-_Fl + dmpcm * _Gl) * _Rp * _Rm * expbeta )
    end
    real(det([M N; P Q]))
end
