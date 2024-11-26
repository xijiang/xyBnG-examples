"""
Results for paper I.

- 1 vs 29 chromosomes
- pg (ulimit -s 1000000) vs tskit
- 5 vs 15 generations
- ΔF = 1.0%, 0.75%, 0.625%, 0.5%
- remove ri-, add ig-, not BLUPs    

Uses paper-i.jl.
"""
function T2024112301(; nrng = 5, nrpt = 50, group = 1)
    clng = [
        1.58534110,
        1.36231102,
        1.21005158,
        1.20000601,
        1.20089316,
        1.17806340,
        1.10682743,
        1.13319770,
        1.05454467,
        1.03308737,
        1.06982474,
        0.87216183,
        0.83472345,
        0.82403003,
        0.85007780,
        0.81013979,
        0.73167244,
        0.65820629,
        0.63449741,
        0.71974595,
        0.69862954,
        0.60773035,
        0.52498615,
        0.62317253,
        0.42350435,
        0.51992305,
        0.45612108,
        0.45940150,
        0.51098607,
    ]
    if group == 1
        paper_1_pg([1.58534110], 0.005, nrng, "rst/pg/c01/1"; nrpt = nrpt)
        paper_1_pg([1.58534110], 0.01,  nrng, "rst/pg/c01/4"; nrpt = nrpt)
    else
        paper_1_pg(clng, 0.005, nrng, "rst/pg/c29/1"; nrpt = nrpt)
        paper_1_pg(clng, 0.01,  nrng, "rst/pg/c29/4"; nrpt = nrpt)
    end
end
