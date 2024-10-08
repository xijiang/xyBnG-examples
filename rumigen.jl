using xyBnG
import xyBnG.Sum: xysum, savesum
import xyBnG.xyTypes: Plan
import xyBnG.xps: gblup, ablup, iblup, aaocs, iiocs, ggocs, agocs, igocs
import xyBnG.xps: savepar, initPop

"""
    rumigen(;
        data = "rst",
        baseDir = "tskit",
        testDir = "cattle",
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
        keep = true,
    )
Function to repeat the simulation for the rumigen project. This is to verify if
my new codes are right. Several issues need to pay attention to.
    
Previously, I forgot to mask the male phenotypes. The number of cadidates, i.e.,
`c` is non-zero, will have more males than females. This is because that EBV for
males are not as accurate as those of females. Male EBV has thus smaller
vairance. 

In the new version of `xyBnG`, I was trying to limit the number of males as set
in `plan`. So I have two `plan`s in this function. If OCS, then `plan` 2 asks
for 50 males and 50 females. The real numbers selected in the simulation are
around 30. `plan` 2 guanruntees the candidates are to be selected.

The other issue is about the GRM. Before I was using allele frequencies accross
generations. Here I only use that of the generation when directional selection
started. Hence selection schemes involving GRM will be different. To verify the
codes, only OCS without GRM are to be considered.
"""
function rumigen(;
    base = "tskit",
    rst = "rumigen/test",
    species = Cattle(5_000),
    trait = Trait("growth", 0.25, 10000),
    nchp = 50_000,
    nref = 10_000,
    nrng = 5,
    nsel = 10,
    plan = Plan(25, 50, 200),
    fixed = ["grt"],
    dF = 0.011,
    maf = 0.0,
    nrpt = 1,
)
    # Scenario recording
    isdir("$rst") && rm("$rst", force = true, recursive = true)
    mkpath("$rst")
    CULLS = (gblup, ablup, iblup)
    OCSS = (aaocs, iiocs, ggocs, agocs, igocs)
    scenario = (
        Data = base,
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
        MAF = maf,
        Nrpt = nrpt,
        Schemes = union(CULLS, OCSS),
    )
    savepar(scenario, "$test/scenario.par")
    chkbase(base, species) # Prepare/verify a base population
    fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"
    pln2 = Plan(50, 50, 200)
    npd = ndigits(nrpt)

    # Simulations
    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        for scheme in CULLS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end
        for scheme in OCSS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(test, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0)
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end
        #=
        if !keep
            for f in readdir("$test")
                occursin(Regex("^$(tag)"), f) && rm("$test/$f", force = true)
            end
        end
        =#
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
