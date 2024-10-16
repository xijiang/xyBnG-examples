using DataFrames
using Serialization
using xyBnG
import xyBnG.Sum: xysum, savesum, cormat, savepar
import xyBnG.xyTypes: Plan
import xyBnG.xps: initPop, chkbase, ggocs, aaocs, iiocs, igocs

"""
Results for paper I.

- 1 vs 29 chromosomes
- pg (ulimit -s 1000000) vs tskit
- 5 vs 15 generations
- ΔF = 1.1% 0.5%
- remove ri-, add ig-, not BLUPs    
"""
function paper_1_tk(bdir; nrpt = 100, ε = 1e-6)
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
    OCSS = (aaocs, ggocs, iiocs, igocs)
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
end

function paper_1_pg(bdir; nrpt = 100, ε = 1e-6)
    # remember `ulimit -s 1000000` before julia REPL
    sdir = @__DIR__
    if !isfile("$sdir/pgsnp")
        @info "Compile pgsnp.cpp"
        run(`g++ -Wall -O3 --std=c++11 -pthread -o $sdir/pgsnp $sdir/pgsnp.cpp`)
    end
end
