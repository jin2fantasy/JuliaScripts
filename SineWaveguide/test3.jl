# constants
const c = 299792458.0
const μ0 = 4e-7pi
const ϵ0 = 1/(μ0*c^2)
# geometry parameters
const a = 0.56e-3
const b = 0.38e-3
const p = 0.23e-3
const h = 0.23e-3
const g = (b-h)/2
const s = p/2
# # space harmonics number in region I
# const n = -3:3
# # standing wave number n region II and III
# const m = 0:4
# number of steps in region II and III
const NN = 500
# wave number in x
const kx = pi/a
#####################
yn(n) = g + h*n/NN
dn(n) = p - 2p/pi*asin((n-1)/NN)
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
# Fl(f, m, n, x) = (lmn_index(f, m, n) * sinh(lmn(f, m, n) * x) +
#     ~lmn_index(f, m, n) * sin(lmn(f, m, n) * x))
# Fpl(f, m, n, x) = (lmn_index(f, m, n) * cosh(lmn(f, m, n) * x) +
#     ~lmn_index(f, m, n) * cos(lmn(f, m, n) * x))
# Gl(f, m, n, x) = (lmn_index(f, m, n) * cosh(lmn(f, m, n) * x) +
#     ~lmn_index(f, m, n) * cos(lmn(f, m, n) * x))
# Gpl(f, m, n, x) = (lmn_index(f, m, n) * sinh(lmn(f, m, n) * x) -
#     ~lmn_index(f, m, n) * sin(lmn(f, m, n) * x))
function Fl(f, m, n, x)
    e = k(f)^2 - kzmn(m, n)^2 - kx^2
    if e < 0
        return sinh(sqrt(-e)*x)
    else
        return sin(sqrt(e)*x)
    end
end
function Fpl(f, m, n, x)
    e = k(f)^2 - kzmn(m, n)^2 - kx^2
    if e < 0
        return cosh(sqrt(-e)*x)
    else
        return cos(sqrt(e)*x)
    end
end
function Gl(f, m, n, x)
    e = k(f)^2 - kzmn(m, n)^2 - kx^2
    if e < 0
        return cosh(sqrt(-e)*x)
    else
        return cos(sqrt(e)*x)
    end
end
function Gpl(f, m, n, x)
    e = k(f)^2 - kzmn(m, n)^2 - kx^2
    if e < 0
        return sinh(sqrt(-e)*x)
    else
        return -sin(sqrt(e)*x)
    end
end
function Fl2(e, x)
    if e < 0
        return sinh(sqrt(-e)*x)
    else
        return sin(sqrt(e)*x)
    end
end
function Fpl2(e, x)
    if e < 0
        return cosh(sqrt(-e)*x)
    else
        return cos(sqrt(e)*x)
    end
end
function Gl2(e, x)
    if e < 0
        return cosh(sqrt(-e)*x)
    else
        return cos(sqrt(e)*x)
    end
end
function Gpl2(e, x)
    if e < 0
        return sinh(sqrt(-e)*x)
    else
        return -sin(sqrt(e)*x)
    end
end

function Rp(b0, n, m)
    _betan = betan(b0, n)
    _mpip = m*pi/p
    if isapprox(abs(_betan), abs(_mpip))
        m == 0 && return convert(Complex, p)
        return convert(Complex, p/2)
    else
        expbeta = cosd(_betan/pi*180p) + im*sind(_betan/pi*180p)
        return (im * _betan * (1 - expbeta * (-1)^m) /
                (_betan^2 - _mpip^2))
    end
end
function Rm(b0, n, m)
    _betan = betan(b0, n)
    _mpip = m*pi/p
    if isapprox(abs(_betan), abs(_mpip))
        m == 0 && return convert(Complex, p)
        return convert(Complex, p/2)
    else
        expbeta = cosd(_betan/pi*180p) - im*sind(_betan/pi*180p)
        return (-im * _betan * (1 - expbeta * (-1)^m) /
                (_betan^2 - _mpip^2))
    end
end

function evalbmpam_dmpcm(f, m)
    bmpam = -Fpl(f, m, NN, h/NN) / Gpl(f, m, NN, h/NN)
    dmpcm = Fpl(f, m, NN, h/NN) / Gpl(f, m, NN, h/NN)
    kf2 = k(f)^2

    for kk in NN:-1:2
        bmpam_index = (bmpam == -1)
        dmpcm_index = (dmpcm == 1)
        yn1 = yn(kk-1)
        dn1 = dn(kk-1)
        dn0 = dn(kk)
        lmk = lmn(f, m, kk)
        lmk1 = lmn(f, m, kk-1)
        e1 = kf2 - kzmn(m, kk)^2 - kx^2
        e2 = kf2 - kzmn(m, kk-1)^2 - kx^2
        Flk = Fl2(e1, .0)
        Fplk = Fpl2(e1, .0)
        Glk = Gl2(e1, .0)
        Gplk = Gpl2(e1, .0)
        Flk1 = Fl2(e2, h/NN)
        Fplk1 = Fpl2(e2, h/NN)
        Glk1 = Gl2(e2, h/NN)
        Gplk1 = Gpl2(e2, h/NN)

        if bmpam_index
            Ymn1n1a = - dn1 / dn0 / lmk
        else
            Ymn1n1a = (dn1/dn0 *
            (Flk + bmpam * Glk) /
            (lmk * (Fplk + bmpam * Gplk))   )
        end
        if dmpcm_index
            Ymn1n1b = dn1 / dn0 / lmk
        else
            Ymn1n1b = (dn1/dn0 *
            (-Flk + dmpcm * Glk) /
            (lmk * (Fplk - dmpcm * Gplk))   )
        end

        bmpam = ( (Ymn1n1a * lmk1 * Fplk1 - Flk1) /
                (Glk1 - Ymn1n1a * lmk1 * Gplk1)     )
        dmpcm = ( (Ymn1n1b * lmk1 * Fplk1 + Flk1) /
                (Glk1 + Ymn1n1b * lmk1 * Gplk1)     )
    end
    (bmpam, dmpcm)::Tuple{Float64, Float64}
end

function detMNPQ(freq, beta0)
    n = -3:3
    np = -3:3 # for clarity
    m = 0:4
    M = Array{Complex{Float64}}(length(n), length(n))
    N, P, Q = similar(M), similar(M), similar(M)
    bmpam_dmpcm = Array(Tuple{Float64, Float64}, length(m))
    for i in 1:length(m)
        bmpam_dmpcm[i] = evalbmpam_dmpcm(freq, m[i])
    end

    for ii in 1:length(np), jj in 1:length(n)
        _nu = nu(freq, beta0, np[ii])
        _Fnu = Fnu(freq, beta0, n[jj])
        _Fpnu = Fpnu(freq, beta0, np[ii])
        _Gnu = Gnu(freq, beta0, n[jj])
        _Gpnu = Gpnu(freq, beta0, np[ii])
        expbeta = exp(im * (betan(beta0, np[ii]) - betan(beta0, n[jj])) * s)

        M[jj,ii] = (ii == jj) ? p^2.0 * _nu * _Fpnu : .0
        N[jj,ii] = (ii == jj) ? p^2.0 * _nu * _Gpnu : .0
        P[jj,ii] = (ii == jj) ? p^2.0 * _nu * _Fpnu : .0
        Q[jj,ii] = (ii == jj) ? -p^2.0 * _nu * _Gpnu : .0
        for kk in 1:length(m)
            _lmn = lmn(freq, m[kk], 1)
            _Fl = Fl(freq, m[kk], 1, .0)
            _Fpl = Fpl(freq, m[kk], 1, .0)
            _Gl = Gl(freq, m[kk], 1, .0)
            _Gpl = Gpl(freq, m[kk], 1, .0)
            _Rp = Rp(beta0, np[ii], m[kk])
            _Rm = Rm(beta0, n[jj], m[kk])

            bmpam, dmpcm = bmpam_dmpcm[kk][1], bmpam_dmpcm[kk][2]
            M[jj,ii] -= (_Fnu * 2 / (1 + (m[kk] == 0)) * _lmn *
            (_Fpl + bmpam * _Gpl) / (_Fl + bmpam * _Gl) * _Rp * _Rm )
            N[jj,ii] -= (_Gnu * 2 / (1 + (m[kk] == 0)) * _lmn *
            (_Fpl + bmpam * _Gpl) / (_Fl + bmpam * _Gl) * _Rp * _Rm )
            P[jj,ii] += (_Fnu * 2 / (1 + (m[kk] == 0)) * _lmn *
            (_Fpl - dmpcm * _Gpl) / (-_Fl + dmpcm * _Gl) * _Rp * _Rm * expbeta )
            Q[jj,ii] -= (_Gnu * 2 / (1 + (m[kk] == 0)) * _lmn *
            (_Fpl - dmpcm * _Gpl) / (-_Fl + dmpcm * _Gl) * _Rp * _Rm * expbeta )
        end
    end
    real(det([M N; P Q]))
end
