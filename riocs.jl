using DataFrames
using LinearAlgebra
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, riocs

function tsbs(;
    nchp = 600,
    nref = 1_100,
    base = "base/chr-1/tskit",
    rst = "1-chr-tskit",
    trait = Trait("growth", 0.25, 1_000),
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
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
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
        Schemes = riocs,
        Nrpt = nrpt,
    )
    savepar(scenario, "$test/scenario.par")
    chkbase(base, species) # Prepare/verify a base population
    npd = ndigits(nrpt)
    
    for irpt = 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        fxy, fmp = "$base/$(species.name).xy", "$base/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        foo, bar = "$tag-rand", tag * '-' * string(riocs)
        riocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end

function pgbs(chr;
    nchp = 600,
    nref = 1_000,
    rst = "1-chr-pgsnp",
    trait = Trait("growth", 0.25, 1_000),
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
    species = Cattle(ne)
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    rst = abspath(rst)

    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        Result = rst,
        Species = species,
        History = hist,
        MutationRate = mr, # per 1e8 bp per meiosis
        Trait = trait,
        Chr_len = chr,
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        MAF = maf,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = riocs,
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
        cmd = `$sdir/pgsnp $(species.nid) $hist $chr $mr`
        run(pipeline(cmd, stdout = "$rst/chr.1", stderr = devnull))
        xyBnG.Conn.PG.toxy("$rst")
        fxy = "$rst/$(species.name).xy"
        fmp = "$rst/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        lmp.chip = lmp.chip .&& .!lmp[!, trait.name] .&& .!lmp[!, "dark"]
        lmp.dark = lmp.dark .&& .!lmp[!, trait.name]

        foo, bar = "$tag-rand", "$tag-riocs"
        riocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum("$rst/summary.ser", summary)
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end

function run_riocs()
    tsbs(nrpt = 20, nchp =  600, rst = "riocs/ts/5d")
    tsbs(nrpt = 20, nchp = 5100, rst = "riocs/ts/5k")
    pgbs(1.585; nrpt = 20, nchp =  600, rst = "riocs/pg/5d")
    pgbs(1.585; nrpt = 20, nchp = 5100, rst = "riocs/pg/5k")
end
