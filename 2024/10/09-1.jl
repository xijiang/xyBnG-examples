"""
Both the random selection and OCS stages uses
`plnb`, which is `Plan(50, 50, 200)`.

## requires
- quick-pgsnp.jl
- xyBnG.jl v1.3.2

## commit
- 
"""
function trial()
    # 2024-10-03
    # includet("quick-pgsnp.jl")
    # commit: 
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

    for i in [1:3; 4:5:29]
        qpg(clng[1:i], rst = "var-length", nrpt = 10)
    end
end