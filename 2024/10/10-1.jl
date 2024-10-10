using Plots

"""
    tstfn(; nsir = 25, ndam = 50, ng = 10, nrpt = 30)
Test ΔF in random mating.

uses "qg/random-mating.jl".
"""
function tstfn(; nsir = 25, ndam = 50, ng = 10, nrpt = 30)
    F = zeros(ng)
    for irpt in 1:nrpt
        irpt % 10 == 0 && print(' ', irpt)
        F += random_mate(nsir, ndam, ng)
    end
    println()
    mF = F / nrpt
    Ne = 2 / mean(x -> 1/x, [nsir, ndam])
    ΔF = 1 / 2Ne
    eF = 1 .- (1 - ΔF) .^ (0:ng-1)
    #round.([mF, ΔF * 0:ng-1], digits = 4)
    [mF eF]
end
