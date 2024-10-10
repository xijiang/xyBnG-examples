using DataFrames
using LinearAlgebra
using Random
using Statistics
using StatsBase

"""
    nrm(ped::DataFrame; m = 30_000)

Given a pedigree `ped`, this function returns a full numerical relationship
matrix, `A`. This function is better for small pedigrees, and for demonstration
only. The maximal matrix size is thus limited to 30k. One can try to set `m` to
a bigger value if your RAM is enough.
"""
function nrm(ped::DataFrame)
    N = size(ped, 1)
    A = zeros(N, N) + I(N)
    for (id, sire, dam) in eachrow(select(ped, :id, :sire, :dam))
        sire * dam ≠ 0 && (A[id, id] += 0.5A[sire, dam])
        for jd = 1:id-1
            sire ≠ 0 && (A[id, jd] = 0.5A[jd, sire])
            dam ≠ 0 && (A[id, jd] += 0.5A[jd, dam])
            A[jd, id] = A[id, jd]
        end
    end
    A
end

"""
    random_mate(nsir::Int, ndam::Int, ng::Int)
Random mate `nsir` sires and `ndam` dams for `ng` generations. Each generation
has `nsir` sires and `ndam` dams. 

This is to confirm ΔF == 1 / 2Ne.
"""
function random_mate(nsir::Int, ndam::Int, ng::Int)
    nid = nsir + ndam
    ped = DataFrame(grt = 0, id = 1:nid, sire = 0, dam = 0)
    ped.sex = shuffle([ones(Int, nsir); zeros(Int, ndam)])
    for ig = 1:ng
        pg = filter(row -> row.grt == ig - 1, ped)
        fra, til = pg.id[end] + 1, pg.id[end] + nid
        sires = pg.id[pg.sex .== 1]
        dams = pg.id[pg.sex .== 0]
        cg = DataFrame(
            grt = ig,
            id = fra:til,
            sire = sample(sires, nid),
            dam = sample(dams, nid),
            sex = shuffle([ones(Int, nsir); zeros(Int, ndam)]),
        )
        append!(ped, cg)
    end
    A = nrm(ped)
    ped.F = diag(A) .- 1
    combine(groupby(ped, :grt), :F => mean => :mF).mF[2:end]
end
