module FSim
function insert!()
end

function splice!()
end

function comb()
end

function ooput()
end

function sim(ne::Int, ng::Int, cl::Float64, mr::Float64)
    # check parameters
    ne < 4 || ng < 100 || cl < 1e-6 || cl > 5 || mr < 0.01 && error("Arguments out of range")

    # Allocate memory
    snp, nsp = zeros(Int8, 120_000, 2, ne, 2), zeros(Int, ne, 2, 2)
    pma, off = 1, 2 # genes flow from 1 to 2
    loci = Dict{Int,Int}()
    λₘ = cl * mr
    pts = begin # maximal mutations per ID
        ept = 2ne * λₘ # expected mutations per generation
        mpt = ept + 10sqrt(ept) # E + 10σ
        zeros(Int, Int(round(mpt)))
    end
    ttbp = Int(1e8 * cl)
    pmt, prc = Poisson(λₘ), Poisson(cl)
    return λₘ
end

end
