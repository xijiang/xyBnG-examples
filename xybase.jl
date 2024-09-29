using DataFrames
using Distributions
using Random

function mutate(nid, rr, mr, bpm, chr...)

end

function drop() end

function xybase(
    nid::Int,
    rr::Float64,
    mr::Float64,
    hist::Int,
    chr::Float64...;
    bpm = 10_000_000,
)
    nid ≥ 4 && rr > 0.0 && mr > 0.0 && hist > 50 && all(chr .> 0.0001) ||
        error("Parameters out of range")
    lmp = DataFrame(chr = Int8[], pos = Int32[])
    tmg = sum(chr) * nid * 2
    tbp = Int(floor(tmg * bpm))
    step = 2000

    for ih ∈ 1:hist
    end
end
