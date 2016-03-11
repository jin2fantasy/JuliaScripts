using Roots
function mysolver(func::Function, lo, hi, N=100)
    guesses = linspace(lo, hi, N)
    f_guess = map(func, guesses)
    a = f_guess[1]
    for i in 2:N
        b = f_guess[i]
        if a * b < 0
            return fzero(func, [guesses[i-1], guesses[i]])
        else
            a = b
        end
    end
    N > 1600 && error("No bracket found/Too many iterations")
    mysolver(func, lo, hi, 2N)
end
