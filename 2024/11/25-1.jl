using DataFrames
using Serialization
using xyBnG

function resum(
    dir::AbstractString,
    sf::AbstractString,
    rpts::UnitRange{Int64},
    ss::AbstractString...;
    pad = 2,
)
    trait = Trait("growth", 0.25, 1000)
    for i in rpts
        tag = lpad(i, pad, '0')
        lmp = deserialize("$dir/$tag-founder.lmp")
        for s in ss
            sm = xyBnG.Sum.xysum("$dir/$tag-$s.ped", "$dir/$tag-$s.xy", lmp, trait)
            xyBnG.Sum.savesum(sf, sm)
        end
        print(" $i")
    end
    println()
end

function merge_repeats(
    til::AbstractString,
    fra::AbstractString,
    rpts::UnitRange{Int64};
    pad = 2,
)
    pr = begin
        file = filter(x -> occursin("aaocs.xy", x), readdir(til))[end]
        parse(Int, split(file, "-")[1])
    end
    fs = readdir(fra)
    for  i in rpts
        tag = lpad(i, pad, '0')
        t2g = lpad(i + pr, pad, '0')
        src = filter(x -> occursin(tag, x), fs)
        for s in src
            tgt = replace(s, tag => t2g)
            mv(joinpath(fra, s), joinpath(til, tgt))
        end
    end
end
