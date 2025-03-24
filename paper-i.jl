using DataFrames
using LinearAlgebra
using Serialization
using Statistics
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, aaocs, iiocs, igocs, hgocs, hhocs, gblup
import xyBnG.xps: F0A, F0G, F0I, F0H

function paper_1_tk(bdir, dF, nrng, rst; nrpt = 100, ε = 1e-6)
    plan = Plan(25, 50, 200; mate = :random)
    plnb = Plan(75, 75, 200; mate = :random)
    fixed = ["grt"]
    nsel, maf = 30, 0.0
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
        fxy = "$bdir/BosTau.xy"
        fmp = "$bdir/BosTau.lmp"
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

function paper_1_pg(chr, dF, nrng, rst; nrpt = 100)
    # remember `ulimit -s 1000000` before julia REPL
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    nchr = length(chr)
    ε = nchr < 5 ? 1e-6 : 0.0
    plan = Plan(25, 50, 200; mate = :random)
    plnb = Plan(75, 75, 200; mate = :random)
    fixed = ["grt"]
    nsel, maf, hist, mr = 30, 0.0, 5000, 4.0
    species = Cattle(plan.noff)
    lgnm = 24.89385779     # length of genome in Morgen
    rog = sum(chr) / lgnm  # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    OCSS = (aaocs, iiocs, ggocs, igocs, hgocs, hhocs)
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
    sumfile = "$rst/summary.ser"
    ocss = (iiocs, ggocs, igocs)
    hocs = (hgocs, hhocs)

    F0 = Dict(  # intialize F0 dictionary
        'a' => 0.0,
        'g' => 0.0,
        'h' => 0.0,
        'i' => 0.0,
    )
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
        lmp, _ = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        begin # update F0
            ped = deserialize("$rst/$tag-rand.ped")
            rxy = "$rst/$tag-rand.xy"
            F0['a'] = F0A(ped)
            F0['g'] = F0G(rxy, lmp, ped)
            F0['i'] = F0I(rxy, lmp, ped)
            F0['h'] = F0H(rxy, lmp, ped)
        end

        # aaocs
        foo, bar = "$tag-rand", "$tag-aaocs"
        aaocs(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0['a'])
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum(sumfile, summary)

        for scheme in ocss
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            ch = string(scheme)[1]
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0[ch]; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
        for scheme in hocs
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0['h']; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end

function pge(chr, dF, nrng, rst; nrpt = 100)
    # remember `ulimit -s 1000000` before julia REPL
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    nchr = length(chr)
    ε = nchr < 5 ? 1e-6 : 0.0
    plan = Plan(25, 50, 200; mate = :random)
    plnb = Plan(75, 75, 200; mate = :random)
    fixed = ["grt"]
    nsel, maf, hist, mr = 30, 0.0, 5000, 4.0
    species = Cattle(plan.noff)
    lgnm = 24.89385779     # length of genome in Morgen
    rog = sum(chr) / lgnm  # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    OCSS = (aaocs, iiocs, ggocs, igocs, hgocs, hhocs)
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
    sumfile = "$rst/summary.ser"
    ocss = (iiocs, ggocs, igocs)
    hocs = (hgocs, hhocs)

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
        ped = deserialize("$rst/$tag-rand.ped")
        rxy = "$rst/$tag-rand.xy"
        Fh = F0H(rxy, lmp, ped)

        for scheme in ocss
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, F0; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
        for scheme in hocs
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plnb, dF, Fh; ε = ε)
            summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            savesum(sumfile, summary)
        end
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end

function cov_g(chr, nrng, rst; nrpt = 100)
    # remember `ulimit -s 1000000` before julia REPL
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
    nchr = length(chr)
    ε = nchr < 5 ? 1e-6 : 0.0
    plan = Plan(25, 50, 200; mate = :random)
    fixed = ["grt"]
    nsel, maf, hist, mr = 30, 0.0, 5000, 4.0
    species = Cattle(plan.noff)
    lgnm = 24.89385779     # length of genome in Morgen
    rog = sum(chr) / lgnm  # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
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
        Nrpt = nrpt,
        Schemes = [gblup],
    )
    savepar(scenario, "$rst/scenario.par")

    npd = ndigits(nrpt)
    open("$rst/desc.txt", "w") do io
        println(io, "BosTau")
        println(io, species.nid)
    end
    sumfile = "$rst/summary.ser"

    for irpt in 1:nrpt
        tag = lpad(irpt, npd, '0')
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        Threads.@threads for i in 1:nchr
            cmd = `$sdir/pgsnp $(species.nid) $hist $(chr[i]) $mr`
            run(pipeline(cmd, stdout = "$rst/chr.$i", stderr = devnull))
        end
        xyBnG.Conn.PG.toxy("$rst")
        fxy = "$rst/$(species.name).xy"
        fmp = "$rst/$(species.name).lmp"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)

        # gblup
        foo, bar = "$tag-rand", "$tag-gblup"
        gblup(rst, foo, bar, lmp, nsel, trait, fixed, plan, ε = ε)
        summary = xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
        savesum(sumfile, summary)
    end
    open("$rst/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
