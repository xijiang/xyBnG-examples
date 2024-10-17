using DataFrames
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, aaocs, iiocs, igocs

function paper_1_tk(bdir, dF, rst; nrpt = 100, ε = 1e-6)
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    nrng, nsel, maf = 5, 30, 0.0
    species = Cattle(5_000)
    lgnm = 24.89385779     # length of genome in Morgen
    blm = deserialize("$bdir/BosTau.lmp")
    lng = combine(groupby(blm, :chr), :pos => maximum => :len)
    rog = sum(lng.len) / lgnm / 1e8 # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    OCSS = (aaocs, ggocs, iiocs, igocs)
    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        BaseDir = bdir,
        Result = rst,
        Species = species,
        Trait = trait,
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        MAF = maf,
        ΔF = dF,
        ε = ε,
        Schemes = OCSS,
        Nrpt = nrpt,
    )
    savepar(scenario, "$rst/scenario.par")
    chkbase(bdir, species) # Prepare/verify a base population
    npd = ndigits(nrpt)
    sumfile = "$rst/summary.ser"
    ocss = (iiocs, ggocs, igocs)
    for irpt ∈ 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        fxy = "$base/BosTau.xy"
        fmp = "$base/BosTau.lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        
        # aaocs
        foo, bar = "$tag-rand", "$tag-aaocs"
        aaocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum(sumfile, summary)

        for scheme in ocss
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end

function paper_1_pg(chr, dF, rst; nrpt = 100)
    # remember `ulimit -s 1000000` before julia REPL
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    nchr = length(chr)
    ε = nchr < 5 ? 1e-6 : 0.0
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    nrng, nsel, maf, hist, mr = 5, 30, 0.0, 5000, 4.0
    species = Cattle(plan.noff)
    lgnm = 24.89385779     # length of genome in Morgen
    rog = sum(chr) / lgnm  # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    OCSS = (aaocs, iiocs, ggocs, igocs)
    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        Result = rst,
        Species = species,
        History = hist,
        MutationRate = mr, # per 1e8 bp per meiosis
        Trait = trait,
        Nchr = nchr,
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
        aaocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum(sumfile, summary)

        for scheme in ocss
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
