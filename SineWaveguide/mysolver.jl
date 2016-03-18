using Roots
function mysolver(func::Function, lo, hi, N=100)
    guesses = linspace(lo, hi, N)
    a = func(guesses[1])
    absminvalue = abs(a)
    absminvalueindex = 1
    for i in 2:N
        b = func(guesses[i])
        if abs(b) < absminvalue
            absminvalue = abs(b)
            absminvalueindex = i
        end
        ratio = b/a
        if ratio < 0
            return fzero(func, [guesses[i-1], guesses[i]])
        else
            a = b
        end
    end
    N > 20000 && error("No bracket found/Too many iterations")
    if N > 6000
        maxindex = min(N, absminvalueindex + 500)
        minindex = max(1, absminvalueindex - 500)
        mysolver(func, guesses[minindex], guesses[maxindex], 2N)
    end
    mysolver(func, lo, hi, 2N)
end
