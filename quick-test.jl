using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, riocs

function qtest(;
    nchp = 600,
    nref = 1_100,
    base = "base/chr-1/tskit",
    rst = "quick-test",
    trait = Trait("growth", 0.25, 1_000),
    plan = Plan(50, 50, 200),
    nrpt = 1,
    nsel = 30,
    nrng = 5,
    maf = 0.0,
    ε = 1e-6,
)
    base, test = abspath(base), abspath(rst)
    isdir(test) && rm(test, force = true, recursive = true)
    mkpath(test)
    species = Cattle(5000)
    fixed = ["grt"]
    dF = 0.011
    OCSS = (riocs, ggocs)
    scenario = (
        BaseDir = base,
        TestDir = test,
        Species = species,
        Trait = trait,
        Nchp = nchp,
        Nref = nref,
        MAF = maf,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plnb,
        Fixed = fixed,
        ΔF = dF,
        ε = ε,
        Schemes = OCSS,
        Nrpt = nrpt,
    )
    savepar(scenario, "$test/scenario.par")
    chkbase(base, species) # Prepare/verify a base population
    npd = ndigits(nrpt)
    sumfile = "$test/summary.ser"

    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, test, plnb, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        # ggocs
        foo, bar = "$tag-rand", "$tag-ggocs"
        ggocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # riocs
        foo, bar = "$tag-rand", "$tag-riocs"
        riocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)
    end

    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
