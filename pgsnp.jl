using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ablup, gblup, iblup, aaocs, ggocs, iiocs, agocs, igocs

"""
    pgsnp(; data = "pgsnp", rst = "rst", nrpt = 1)
# Using `pgsnp` results as base population.

## Manule

Suppose you are in a working directory ``dir``. This file ``pgsnp.jl`` is in
``script/``. You want to save the results in ``rst``

- Date: 2024-08-27, works ok with xyBnG v1.2.4
"""
function pgsnp(;
    nchp = 600,
    nref = 1_000,
    rst = "1-chr-pgsnp",
    trait = Trait("growth", 0.25, 1_000),
    nchr = 1,
    nrpt = 1,
    nsel = 30,
    nrng = 5,
    ne = 200,
    maf = 0.0,
    dF = 0.011,
    hist = 5000,
    mr = 4.0,
)
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    1 ≤ nchr ≤ 29 || (nchr = 1)
    species = Cattle(ne)
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    CULLS = (gblup, ablup, iblup)
    OCSS = (aaocs, iiocs, ggocs, agocs, igocs)
    rst = abspath(rst)

    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        Result = rst,
        Species = species,
        History = hist,
        MutationRate = mr, # per 1e8 bp per meiosis
        Trait = trait,
        Nchr = nchr, # = 1:29
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        MAF = maf,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = union(CULLS, OCSS),
    )
    savepar(scenario, "$rst/scenario.par")

    npd = ndigits(nrpt)
    chrln = [
        1.58,
        1.36,
        1.21,
        1.20,
        1.18,
        1.11,
        1.13,
        1.05,
        1.03,
        1.07,
        0.87,
        0.83,
        0.82,
        0.85,
        0.81,
        0.73,
        0.66,
        0.63,
        0.71,
        0.69,
        0.61,
        0.52,
        0.62,
        0.42,
        0.52,
        0.46,
        0.46,
        0.51,
    ]
    open("$rst/desc.txt", "w") do io
        println(io, "BosTau")
        println(io, species.nid)
    end
    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        Threads.@threads for i ∈ 1:nchr
            cmd = `$sdir/pgsnp $(species.nid) $hist $(chrln[i]) $mr`
            run(pipeline(cmd, stdout = "$rst/chr.$i", stderr = devnull))
        end
        xyBnG.Conn.PG.toxy("$rst")
        fxy = "$rst/$(species.name).xy"
        fmp = "$rst/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        for scheme in CULLS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plan)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum("$rst/summary.ser", summary)
        end
        for scheme in OCSS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum("$rst/summary.ser", summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
