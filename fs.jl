module FDR
using Distributions
using Random

"""
    sim(nid, ngn, cln, mr; rr = 1.0, mm=120_000, bpm=100_000_000, cpd=1_000)
This is forward simulator of one chromosome of a diploid Fisher-Wright
population:

- Populaiton size: `nid`. 
- Chromosome length: `cln` Morgen
- Mutation rate: `mr` per Morgen.

Optional arguments:
- Recombination rate: `rr` = 1.0 crossover per Morgan per meiosis.
- Maximal number of mutations per haploid: `mm` = 1.2e5
- Number of base pairs per Morgen: `bpm` = 1e8
- Comb period: `cpd` = 1e3 generations, to remove fixed loci every `cpd`
  generations.
"""
function sim(
    nid::Int,       # population size
    ngn::Int,       # number of generations
    cln::Float64,   # chromosome length in Morgen
    mr::Float64;    # mutation rate per Morgen per generation
    rr = 1.0,       # recombination rate per Morgen per meiosis
    mm = 120_000,   # maximal mutations per haploid
    bpm = 100_000_000,  # Base pair per Morgen
    comb = 1_000,    # comb (out fixed loci) period in generations
)
    # check parameters
    nid < 4 ||
        ngn < 100 ||
        cln < 1e-6 ||
        cln > 5 ||
        mr < 0.01 && error("Arguments out of range")

    @info "Generating an ideal population"
    @info "  - Initializing pars and storage"
    # Allocate memory
    snp, nsp = zeros(Int8, mm, 2, nid, 2), zeros(Int, nid, 2, 2)
    pma, off = 1, 2 # genes flow from 1 to 2
    loci = Dict{Int,Int}()
    λₘ = cln * mr
    pts = begin # maximal mutations per ID
        ept = 2nid * λₘ # expected mutations per generation
        mpt = ept + 10sqrt(ept) # E + 10σ
        zeros(Int, Int(round(mpt)))
    end
    ttbp = Int(bpm * cln)
    pmt, prc = Poisson(λₘ), Poisson(cln)
    return λₘ
    @info "  - Balancing mutations and recombinations ..."
    for ign ∈ 1:ngn
        for ihp ∈ 1:2
            nmt = rand(pmt)
            mts = rand(1:ttbp, nmt)
        end
    end
end
end # module FS