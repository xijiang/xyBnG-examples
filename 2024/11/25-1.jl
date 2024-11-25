using DataFrames
using Serialization
using xyBnG

function resum(dir::AbstractString)
    ss = ("aaocs", "ggocs", "igocs", "iiocs")
    trait = Trait("growth", 0.25, 1000)
    for i in 1:50
        tag = lpad(i, 2, '0')
        lmp = deserialize("$dir/$tag-founder.lmp")
        for s in ss
            sm = xyBnG.Sum.xysum("$dir/$tag-$s.ped", "$dir/$tag-$s.xy", lmp, trait)
            xyBnG.Sum.savesum("summary.ser", sm)
        end
        print(" $i")
    end
    println()
end
