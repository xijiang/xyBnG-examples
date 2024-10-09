using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, iiocs, aaocs, riocs

"""
    qpg(; data = "pgsnp", rst = "rst", nrpt = 1)
# Using `pgsnp` results as base population.

## Manual

Suppose you are in a working directory ``dir``. This file ``quick-pgsnp.jl`` is
in ``dir/script/``. You want to save the results in ``dir/rst/``

Notes: the chip, reference and QTL SNPs are sampled independently, meaning they
might share a few common loci.

The number of QTL and reference loci are 401.7/Morgen. This is assume a cattle
has a genome of 24.894 M. We assume there are 50k chip SNP, 10k SNP for
reference and QTL in the genome. The heritability is adjusted correspondingly.
"""
function qpg(chr; rst = "", nrpt = 1)
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    # Simulation parameters
    nchr = length(chr)
    ε = nchr == 1 ? 1e-6 : 0.0
    plan = Plan(25, 50, 200)
    nrng, nsel, maf, dF, hist, mr = 5, 30, 0.0, 0.011, 5000, 4.0
    species = Cattle(plan.noff)
    lgnm = 24.89385779     # length of genome in Morgen
    rog = sum(chr) / lgnm  # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    plnb, fixed = Plan(50, 50, 200), ["grt"]  # plnb is in case
    OCSS = (aaocs, iiocs, ggocs, riocs)

    rst = joinpath(abspath(rst), "$nchr")
    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)

    scenario = (
        Result = rst,
        Species = species,
        History = hist,
        MutationRate = mr, # per 1e8 bp per meiosis
        Trait = trait,
        Nchr = length(chr),
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        MAF = maf,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = OCSS,
    )
    savepar(scenario, "$rst/scenario.par")

    npd = ndigits(nrpt)
    open("$rst/desc.txt", "w") do io
        println(io, "BosTau")
        println(io, species.nid)
    end
    
    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        Threads.@threads for i = 1:nchr
            cmd = `$sdir/pgsnp $(species.nid) $hist $(chr[i]) $mr`
            run(pipeline(cmd, stdout = "$rst/chr.$i", stderr = devnull))
        end
        xyBnG.Conn.PG.toxy("$rst")
        fxy = "$rst/$(species.name).xy"
        fmp = "$rst/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)

        # aaocs
        foo, bar = "$tag-rand", "$tag-aaocs"
        ggocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # ggocs
        foo, bar = "$tag-rand", "$tag-ggocs"
        ggocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # riocs
        foo, bar = "$tag-rand", "$tag-riocs"
        riocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)

        # iiocs
        foo, bar = "$tag-rand", "$tag-iiocs"
        riocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
