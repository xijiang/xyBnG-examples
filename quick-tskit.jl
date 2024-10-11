using DataFrames
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, riocs, aaocs, iiocs


function qtskit(cdir::Int; nrpt = 30)
    ε = cdir < 5 ? 1e-6 : 0.0
    plan, plnb, fixed = Plan(25, 50, 200), Plan(50, 50, 200), ["grt"]
    nrng, nsel, maf, dF = 5, 30, 0.0, 0.011
    species = Cattle(5_000)
    lgnm = 24.89385779     # length of genome in Morgen
    base = "base/vchr/$cdir"
    blm = deserialize("$base/BosTau.lmp")
    lng = combine(groupby(blm, :chr), :pos => maximum => :len)
    rog = sum(lng.len) / lgnm / 1e8 # ratio of genome
    nchp = Int(round(5e4 * rog))
    nref = Int(round(1e4 * rog))
    trait = Trait("growth", 0.25 * rog, nref) # nqtl = nref
    OCSS = (aaocs, iiocs, ggocs, riocs)
    rst = "rst/vchr/$cdir"
    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    scenario = (
        BaseDir = base,
        Result = rst,
        Species = species,
        Trait = trait,
        Nchr = cdir,
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
    chkbase(base, species) # Prepare/verify a base population
    npd = ndigits(nrpt)
    sumfile = "$rst/summary.ser"
    OCSS = (iiocs, ggocs, riocs)
    
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

        for scheme in OCSS
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
