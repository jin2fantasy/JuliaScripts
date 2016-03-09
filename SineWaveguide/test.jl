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
freq = 190e9
beta0 = 0.0
#####################
yn(n) = g + h*n/NN
dn(n) = p - p/pi*acos((NN-2n+2)/NN)
k = 2pi*freq/c
kzmn(m, n) = m*pi / dn(n)
betan(n) = beta0 + 2pi*n/p
nu(n) = sqrt(abs(k^2 - betan(n)^2 - kx^2))
nu_index(n) = (k^2 - betan(n)^2 - kx^2) < 0
lmn(m, n) = sqrt(abs(k^2 - kzmn(m ,n)^2 - kx^2))
lmn_index(m, n) = (k^2 - kzmn(m, n)^2 - kx^2) < 0
Fnu(n) = nu_index(n) * sinh(nu(n) * g) + ~nu_index(n) * sin(nu(n) * g)
Fpnu(n) = nu_index(n) * cosh(nu(n) * g) + ~nu_index(n) * cos(nu(n) * g)
Gnu(n) = nu_index(n) * cosh(nu(n) * g) + ~nu_index(n)* cos(nu(n) * g)
Gpnu(n) = nu_index(n) * sinh(nu(n) * g) - ~nu_index(n) * sin(nu(n) * g)
Fl(m, n, x) = (lmn_index(m, n) * sinh(lmn(m, n) * x) +
    ~lmn_index(m, n) * sin(lmn(m, n) * x))
Fpl(m, n, x) = (lmn_index(m, n) * cosh(lmn(m,n) * x) +
    ~lmn_index(m, n) * cos(lmn(m, n) * x))
Gl(m, n, x) = (lmn_index(m, n) * cosh(lmn(m, n) * x) +
    ~lmn_index(m, n) * cos(lmn(m, n) * x))
Gpl(m, n, x) = (lmn_index(m, n) * sinh(lmn(m, n) * x) -
    ~lmn_index(m, n) * sin(lmn(m, n) * x))

function Rp(n, m)
    if isapprox(abs(betan(n)), abs(m*pi/p))
        m == 0 && return p
        return p/2
    else
        return (im * betan(n) * (1 - exp(im * betan(n) * p) * (-1)^m) /
                (betan(n)^2 - (m*pi/p)^2))
    end
end
function Rm(n, m)
    if isapprox(abs(betan(n)), abs(m*pi/p))
        m == 0 && return p
        return p/2
    else
        return (-im * betan(n) * (1 - exp(-im * betan(n) * p) * (-1)^m) /
                (betan(n)^2 - (m*pi/p)^2))
    end
end

function evalbmpam_dmpcm(m)
    bmpam = -Fpl(m, NN, g+h) / Gpl(m, NN, g+h)
    dmpcm = Fpl(m, NN, g+h) / Gpl(m, NN, g+h)
    @inbounds for kk in NN:-1:2
        bmpam_index = (bmpam == -1)
        dmpcm_index = (dmpcm == -1)
        if bmpam_index
            Ymn1n1a = - dn(kk-1) / dn(kk) / lmn(m, kk)
        else
            Ymn1n1a = (dn(kk-1)/dn(kk) *
            (Fl(m, kk, yn(kk-1)) + bmpam * Gl(m, kk, yn(kk-1))) /
            (lmn(m, kk) * (Fpl(m, kk, yn(kk-1)) + bmpam * Gpl(m, kk, yn(kk-1))))   )
        end
        if dmpcm_index
            Ymn1n1b = - dn(kk-1) / dn(kk) / lmn(m, kk)
        else
            Ymn1n1b = (dn(kk-1)/dn(kk) *
            (Fl(m, kk, yn(kk-1)) + dmpcm * Gl(m, kk, yn(kk-1))) /
            (lmn(m, kk) * (Fpl(m, kk, yn(kk-1)) + dmpcm * Gpl(m, kk, yn(kk-1))))   )
        end

        bmpam = ( (Ymn1n1a .* lmn(m, kk-1) .* Fpl(m, kk-1, yn(kk-1)) - Fl(m, kk-1, yn(kk-1))) ./
                (Gl(m, kk-1, yn(kk-1)) - Ymn1n1a .* lmn(m, kk-1) .* Gpl(m, kk-1, yn(kk-1)))     )
        dmpcm = ( (Ymn1n1b .* lmn(m, kk-1) .* Fpl(m, kk-1, yn(kk-1)) - Fl(m, kk-1, yn(kk-1))) ./
                (Gl(m, kk-1, yn(kk-1)) - Ymn1n1b .* lmn(m, kk-1) .* Gpl(m, kk-1, yn(kk-1)))     )
    end
    bmpam, dmpcm
end


n = -3:3
np = -3:3 # for clarity
m = 0:4
M = zeros(length(n), length(n))
N, P, Q = similar(M), similar(M), similar(M)

@inbounds for ii in 1:length(np), jj in 1:length(n), kk in 1:length(m)
    bmpam, dmpcm = evalbmpam_dmpcm(m[kk])
    M[jj,ii] = ( (ii == jj) * p^2 * nu(np[ii]) * Fpnu(np[ii]) -
    Fnu(n[jj]) * 2 / (1 + (kk == 0)) * lmn(kk, 1) * (Fpl(kk, 1, g) +
    bmpam * Gpl(kk, 1, g)) / (Fl(kk, 1, g) + bmpam * Gl(kk, 1, g)) *
    Rp(np[ii], kk) * Rm(n[jj], kk) )
    N[jj,ii] = ( (ii == jj) * p^2 * nu(np[ii]) * Gpnu(np[ii]) -
    Gnu(n[jj]) * 2 / (1 + (kk == 0)) * lmn(kk, 1) * (Fpl(kk, 1, g) +
    bmpam * Gpl(kk, 1, g)) / (Fl(kk, 1, g) + bmpam * Gl(kk, 1, g)) *
    Rp(np[ii], kk) * Rm(n[jj], kk) )
    P[jj,ii] = ( (ii == jj) * p^2 * nu(np[ii]) * Fpnu(np[ii]) +
    Fnu(n[jj]) * 2 / (1 + (kk == 0)) * lmn(kk, 1) * (Fpl(kk, 1, g) -
    dmpcm * Gpl(kk, 1, g)) / (-Fl(kk, 1, g) + dmpcm * Gl(kk, 1, g)) *
    Rp(np[ii], kk) * Rm(n[jj], kk) * exp(im * (betan(np[ii]) - betan(n[jj]) * s)) )
    Q[jj,ii] = ( -(ii == jj) * p^2 * nu(np[ii]) * Gpnu(np[ii]) -
    Gnu(n[jj]) * 2 / (1 + (kk == 0)) * lmn(kk, 1) * (Fpl(kk, 1, g) -
    dmpcm * Gpl(kk, 1, g)) / (-Fl(kk, 1, g) + dmpcm * Gl(kk, 1, g)) *
    Rp(np[ii], kk) * Rm(n[jj], kk) * exp(im * (betan(np[ii]) - betan(n[jj]) * s)) )
end
det([M N; P Q])
