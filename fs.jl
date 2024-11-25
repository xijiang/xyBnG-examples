using Distributions
using Random

"""
    fwdsim(
        nid::Int;       # population size
        ngn::Int = 5000,    # number of generations
        cln::Float64 = 1.0, # chromosome length in Morgen
        mr::Float64 = 1.0,  # mutation rate per Morgen per generation
        rr::Float64 = 1.0,  # recombination rate per Morgen per meiosis
        )
This is forward simulator of one chromosome of a diploid Fisher-Wright
population.
"""
function fwdsim(
    nid::Int;       # population size
    ngn::Int = 5_000,   # number of generations
    cln::Float64 = 1.0, # chromosome length in Morgen
    mr::Float64 = 1.0,  # mutation rate per Morgen per generation
    rr::Float64 = 1.0,  # recombination rate per Morgen per meiosis
    )
    mm = 120_000      # maximal mutations per haploid
    bpm = 100_000_000 # Base pairs per Morgen
    comb = 1_000      # comb (out fixed loci) period in generations
    
    @info "Generating an ideal population"
    @info "Check parameters"
    nid ≥ 4 && ngn ≥ 100 && 1e-6 ≤ cln ≤ 6 && 01 ≤ mr ≤ 12 && 0.1 ≤ rr ≤ 5||
        error("Arguments out of range")

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
