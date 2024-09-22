using DataFrames
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase

function fewer_chrs(;
    nchp = 32_000,
    nref = 1_100,
    base = "base/chr-1/tskit",
    rst = "rst",
    trait = Trait("growth", 0.25, 1_000),
    nrpt = 1,
    nsel = 20,
    nrng = 5,
    maf = 0.0,
    ε=1e-6
)
    script = @__DIR__
    base, test = "$script/../$base", "$script/../$rst"
    isdir(test) && rm(test, force = true, recursive = true)
    mkpath(test)
    species = Cattle(5000)
    plan, fixed = Plan(25, 50, 200), ["grt"]
    dF = 0.011
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
        Plan = plan,
        Fixed = fixed,
        ΔF = dF,
        ε = ε,
        Nrpt = nrpt,
    )
    savepar(scenario, "$test/scenario.par")
    chkbase(base, species) # Prepare/verify a base population
    npd = ndigits(nrpt)
    pln2 = Plan(50, 50, 200)
    sumfile = "$test/summary.ser"

    for irpt in 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        # GBLUP
        foo, bar = "$tag-rand", "$tag-gblup"
        gblup(test, foo, bar, lmp, nsel, trait, fixed, plan; ε = ε)
        summary = xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
        savesum(sumfile, summary)

        # GG-OCS
        foo, bar = "$tag-rand", "$tag-ggocs"
        ggocs(test, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = ε)
        summary = xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
        savesum(sumfile, summary)
        
        #continue
        ## II-OCS
        foo, bar = "$tag-rand", "$tag-iiocs"
        iiocs(test, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0; ε = ε)
        summary = xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
        savesum(sumfile, summary)
    end
end
