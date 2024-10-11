using DataFrames
using Serialization
using xyBnG.XY

function subtskit(; idir = "base/tskit", odir = "base/vchr")
    lmp = deserialize("$idir/BosTau.lmp")
    _, nhp = XY.dim("$idir/BosTau.xy")
    for nchr in [1:3; 4:5:29]
        @info "Processing chromosome $nchr"
        tdir = joinpath(odir, "$nchr")
        isdir(tdir) || mkpath(tdir)
        slc = sum(lmp.chr .≤ nchr)
        cp(joinpath(idir, "desc.txt"), joinpath(tdir, "desc.txt"), force = true)
        serialize(joinpath(tdir, "BosTau.lmp"), lmp[1:slc, :])
        XY.sub("$idir/BosTau.xy", 1:slc, 1:nhp, joinpath(tdir, "BosTau.xy"))
    end
end