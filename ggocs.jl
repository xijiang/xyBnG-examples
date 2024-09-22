"""
    ggocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F0; ε = 1e-6)
Optimal contribution selection with `G` relationship matrix for both EBV and
constraint on `foo`.xy and `foo`.ped in directory `test` for `ngn` generations.
SNP linkage information are in DataFrame `lmp`. The results are saved in
`bar`.xy, `bar`.ped in directory `test`. The selection is on a single trait
`trait` with fixed effects `fixed`, which is a column name vector in pedigree
DataFrame. Parents are sampled according to `plan`. The constraint ΔF is `dF`.
`F0` is the inbreeding coefficient of the `foo` population.

This function uses the TM1997 algorithm for OCS. As the constraint matrix is
using GRM, which has already considered the allele frequency changes, hence
option `ong` is set to `true`.

See also [`randbrd`](@ref), [`aaocs`](@ref), [`iiocs`](@ref), [`iiocs`](@ref),
[`agocs`](@ref), [`igocs`](@ref).
"""
function ggocs(test, foo, bar, lmp, ngn, trait, fixed, plan, dF, F; ε = 1e-6)
    @info "  - Directional selection GGOCS for $ngn generations"
    ped, xy = deserialize("$test/$foo.ped"), "$test/$bar.xy"
    cp("$test/$foo.xy", xy, force=true)
    for ign = 1:ngn
        print(" $ign")
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, trait)
        G = grm(xy, lmp.chip, lmp.frq) + I * ε
        giv = inv(G)
        Predict!(ids, ped, fixed, giv, trait)
        g22 = G[ids, ids]
        ng = Select(ids, plan, ped, g22, trait, dF, ign; F0=F0)
        reproduce!(ng, ped, xy, lmp, trait)
    end
    println()
    serialize("$test/$bar.ped", ped)
end
