using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ablup, gblup, iblup, aaocs, ggocs, iiocs, agocs, igocs

function eaap(;
    rst = "eaap",
    base = "tskit",
    species = Cattle(5_000),
    trait = Trait("growth", 0.25, 10_000),
    nchp = 50_000,
    nref = 10_000,
    nrng = 5,
    nsel = 20,
    plan = Plan(25, 50, 200),
    fixed = ["grt"],
    dF = 0.011,
    nrpt = 1,
) 
    isdir(rst) && rm(rst, force = true, recursive = true)
    base = abspath("../base/$base")
    mkpath(rst)
    CULLS = (gblup, ablup, iblup)
    OCSS = (aaocs, iiocs, ggocs, agocs, igocs)
    scenario = (
        Base = base,
        Results = rst,
        Species = species,
        Trait = trait,
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = union(CULLS, OCSS),
    )
    savepar(scenario, "$rst/scenario.par")
    chkbase(base, species) # Prepare/verify a base population
    fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"
    pln2 = Plan(50, 50, 200)
    npd = ndigits(nrpt)

    for irpt ∈ 1:nrpt
        tag = lpad(irpt, npd, '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        # initPop: sample founder, then random selection for nrng generations
        lmp, F0 = initPop(fxy, fmp, rst, plan, 0.0, nchp, nref, nrng, trait, tag)

        # ablup
        foo, bar = "$tag-rand", "$tag-ablup"
        ablup(rst, foo, bar, lmp, nsel, trait, fixed, plan)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # gblup
        foo, bar = "$tag-rand", "$tag-gblup"
        gblup(rst, foo, bar, lmp, nsel, trait, fixed, plan; ε = 0.)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # iblup
        foo, bar = "$tag-rand", "$tag-iblup"
        iblup(rst, foo, bar, lmp, nsel, trait, fixed, plan; ε=0.0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # aaocs
        foo, bar = "$tag-rand", "$tag-aaocs"
        aaocs(rst, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # ggocs
        foo, bar = "$tag-rand", "$tag-ggocs"
        ggocs(rst, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = 0.0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # iiocs
        foo, bar = "$tag-rand", "$tag-iiocs"
        iiocs(rst, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = 0.0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # agocs
        foo, bar = "$tag-rand", "$tag-agocs"
        agocs(rst, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = 0.0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # igocs
        foo, bar = "$tag-rand", "$tag-igocs"
        igocs(rst, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = 0.0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)
    end
    open("$rst/scenario.par", "w") do io
        println(io, "Ended: ", time())
    end
end
