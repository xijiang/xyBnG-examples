"""
    function ss_brd_cycle(xy, lmp, ped, idxwt, trts...)
The breeding cycle for Sustainable sheep.
"""
function ss_brd_cycle(xy, lmp, ped, idxwt, trts...)
    @info "   .1. Calculate EBV and selection index"
    nhp = size(xy, 2)
    M = view(xy, lmp.chip, 1:2:nhp) + view(xy, lmp.chip, 2:2:nhp)
    G = RS.grm(M; p = lmp.frq[lmp.chip]) + 1e-4I
    giv = inv(G)
    ID = ped.id[ped.grt .== ped.grt[end]]
    Predict!(ID, ped, ["grt"], giv, trts...)

    for t in trts
        ped[ID, :idx] += idxwt[t.name] * ped[ID, "ebv_" * t.name]
    end

    @info "   .2. Select parents"
    @info "       -------- Waiting for rules --------"
    @info "   .3. Generate parent pairs"
    @info "       -------- Waiting for rules --------"
    @info "   .4. Create offspring with trait `litter_size`"
    @info "   .5. Decide AI, elite and candidate rams"
    @info "       -------- Waiting for rules --------"
    @info "   .6. Cull some ewes"
    @info "       -------- Waiting for rules --------"
    @info "   .7. Cull some rams"
    @info "       -------- Waiting for rules --------"
end

"""
    function ssheep(dir::AbstractString)
This is for the digital twin of the NSG sheep population. Lines end with `<==`
comments are parameters to be decided.

The script below is self-illustrative. Please see the comments for details.
"""
function ssheep(dir::AbstractString; new_base = false, nrpt = 1)
    @info "1. Creating a sheep base population"
    nws = Sheep(1000) # <== Norwegian White Sheep
    ngh = 5000        # <== number of generations in history
    mr = 2.0          # <== number of mutations per M per meiosis
    rr = 1.0          # <== number of recombination per M per meiosis
    M = 1e8           # <== base pair per Morgan
    fxy, fmp = joinpath.(dir, ("nws.xy", "nws.lmp"))
    if !all(isfile.([fxy, fmp])) || new_base
        muts, cbp = Fdr.FWP(nws.nid, ngh, nws.chromosome, mr; rr = rr, M = M)
        xy, lmp = Fdr.muts2xy(muts, cbp; flip = true)
        Conn.Mat.toxy(fxy, xy) # save in dir
        serialize(fmp, lmp)
    end # 45 minutes with AMD 3900x, ~7G DDR4-2133 RAM used => ~1.5m loci

    ped = nothing
    @info "2. The breeding cycles"
    for irpt ∈ 1:nrpt
        println()
        @info "\n\tRepeat $irpt / $nrpt"
        println()

        @info "2.1. Sampling ID from the base population"
        ne = 200  # <== effective population size
        maf = 0.0 # (also) remove fixed loci
        xy, lmp = Conn.xy.sampleID(fxy, fmp, maf, ne)

        @info "2.2. Expanding the population exponentially to prod. pop. size"
        ppsz = 500 # <== production population size
        ngrn = 10  # <== number of generations to expand to ppsz
        xy, lmp = expand(xy, lmp, ppsz, ngrn)

        @info "2.3. Sampling chip, ref and QTL loci"
        chip = SNPSet("chip", 50_000; maf = 0.3, exclusive = true) # <==
        dark = SNPSet("dark", 10_000; maf = 0.0, exclusive = true) # <== ref was for allele names
        qtl = SNPSet("qtl", 10_000; maf = 0.05, exclusive = true)  # <==
        xy, lmp = Conn.xy.sampleLoci(xy, lmp, chip, dark, qtl)

        @info "2.4. Creating traits"
        growth = Trait("grθ"; h² = 0.5)     # <==
        carcass = Trait("kks"; h² = 0.3)    # <==
        maternal = Trait("mtn"; h² = 0.1)   # <==
        litter_size = Trait("lsz", [10.0, 20.0, 30.0]; h² = 0.2)  # <==
        trts = [growth, carcass, maternal, litter_size]
        R = [1 0.7 0.01 0.2; 0.7 1 0.05 0.3; 0.01 0.05 1 0.15; 0.2 0.3 0.15 1] # <==

        ped = trait(xy, lmp, R, growth, carcass, maternal, litter_size)
        ped.idx .= 0.0
        # Above creates a population to start with, this takes < 1 minute

        nbc = 1 # <== number of breeding cycles, or, generations
        idxwt = Dict(zip(getfield.(trts, :name), [1.0, 1.0, 1.0, 1.0])) # <== weights of trait EBV in selection index
        for ibc ∈ 1:nbc
            @info "2.5. The breeding cycles: $ibc / $nbc"
            ss_brd_cycle(xy, lmp, ped, idxwt, trts...)
        end
    end

    return xy, lmp, ped
end
