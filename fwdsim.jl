using BenchmarkTools
using DataFrames
using Distributions
using Random

"""
    rndlmp(nlc::Int, cbp::Int...) -> DataFrame
Sample random SNP loci within the given chromosome lengths `cbp...` in base
pairs. The total number of loci is `nlc`. The number of loci on each
chromosome is proportional to its length.
"""
function rndlmp(nlc::Int, cbp::Int...)
    nch = length(cbp)
    bp = ['A', 'C', 'G', 'T']
    tbp = sum(cbp)          # genome length in base pairs
    csc = Int.(round.(cbp ./ tbp .* nlc)) # chromosome snp counts
    sum(csc .<= cbp) == nch ||
        throw(ArgumentError("No. loci must be ≪ the total length in bp"))
    lmp = DataFrame(chr = Int8[], pos = Int32[], ref = Char[], alt = Char[])
    for (chr, (len, cnt)) in enumerate(zip(cbp, csc))
        append!(
            lmp,
            DataFrame(
                chr = chr,
                pos = sort(sample(1:len, cnt, replace = false)),
                ref = 'A',
                alt = 'A',
            ),
        )
    end
    for r in eachrow(lmp)
        r.ref, r.alt = sample(bp, 2, replace = false)
    end
    lmp
end

function anneal(
    nid::S,
    nlc::S,
    ngt::S,
    mtr::T,
    rcr::T,
    cln::T...;
    bpm = 1e8,
) where {S<:Int,T<:Float64}
    cbp = (Int ∘ floor).(cln .* bpm)
    lmp = rndlmp(nlc, cbp...)
    pop = BitArray(undef, nlc, nid, 2)
    rand!(pop)  # pop now has nlc loci with a roughly flat freq spectrum
    alt = BitArray(undef, nlc, nid, 2)
    emt = mtr / bpm * nlc * 2 * nid # expected no. of mutations in `pop`/generation
    nmut = 0
    for igt ∈ 1:ngt
        nmut += rand(Poisson(emt))
        println(nmut)
    end
end

function fwdsim(nlc, nid, ngt, cbp, mtr, rcr; bpm = 1e8)
    nhp, tpt = 2nid, 2nid * nlc
    pop = BitArray(undef, nlc, nhp)
    alt = BitArray(undef, nlc, ngt)
    rand!(pop)  # pop now has nlc loci with a roughly flat freq spectrum
    emt = mtr / bpm * tpt # expected number of mutations in pop
    erc = rcr * cbp
    nrc = zeros(Int, nhp) # number of recombinations happens in each haplotype of alt
    sex = [zeros(nid); ones(nid)]
    for igt = 1:ngt
        # number of mutations
        nmt = rand(Poisson(emt))
        # where are they?
        idx = rand(1:tpt, nmt)
        # flip the bits
        pop[idx] = .!pop[idx]
        # how many recombinations?
        rand!(Poisson(erc), nrc)
        # sexes in pop
        shuffle!(sex)
        # sires in pop
        sires = findall(x -> x == 1, sex)
        # dams in pop
        dams = findall(x -> x == 0, sex)
        # sample parents
        pm = sortslices([sample(sires, nid) sample(dams, nid)], dims = 1)
        # how many recombinations in each haplotype of alt
        rand!(Poisson(erc), nrc)
        # where are they?

    end
end
