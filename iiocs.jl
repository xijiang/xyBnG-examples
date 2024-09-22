"""
    iiocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
Optimal contribution selection with `IBD` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parents are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS.
"""
function iiocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0)
    @info "  - Directional selection IIOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    G = nothing
    if isfile("$test/$foo.irm")
        G = zeros(size(ped, 1), size(ped, 1))
        read!("$test/$foo.irm", G)
    else
        @info "  - Calculating IBD relationship matrix"
        G = irm(xy, lmp.chip, 1:size(ped, 1))
        write("$test/$foo.irm", G)
    end
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        mid = size(ped, 1)
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
        G = xirm(G, xy, lmp.chip, mid, size(ped, 1))
    end
    println()
    serialize("$test/$bar.ped", ped)
end
