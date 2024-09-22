"""
    gblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
Directional selection with `G` relationship matrix for EBV on `foo`.xy and
`foo`.ped in directory `test` for `ngn` generations. SNP linkage information are
in DataFrame `lmp`. The results are saved in `bar`.xy, `bar`.ped in directory
`test`. The selection is on a single trait `trait` with fixed effects `fixed`,
which is a column name vector in pedigree DataFrame. The selection is based on
the selection plan `plan`.

See also [`ablup`](@ref), [`iblup`](@ref).
"""
function gblup(test, foo, bar, lmp, ngn, trait, fixed, plan)
    @info "  - Directional selection GBLUP for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq)
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        ng = Select(ids, plan, ped, trait)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end
