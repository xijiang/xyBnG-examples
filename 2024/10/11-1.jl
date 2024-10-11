"""
using tskit result to test about various chromosome numbers.
"""
function trial()
    # 2024-10-11
    # includet("quick-tskit.jl")
    # commit: 
    for nchr in [1:3; 4:5:29]
        qtskit(nchr; nrpt = 20)
    end
end
