# includet("script/summary.jl")
function gatherdata(dir)
    pg = []
    for i in 1:4
        dat = getdata(joinpath(dir, "$i"))
        push!(pg, dat)
    end
    pg
end

function figprt(; bs = "pg", nrs = "05", chr = "c01")
    pg = gatherdata("/home/xijiang/workspace/xyBnG/rst/paper-1/$bs/$nrs/$chr")
    fdr = -parse(Int, nrs)

    ss = keys(pg[1])
    low, hgh = lmtof(pg, ss, :nsire, :ndam; fg = 2 - fdr)
    ylm = (low -1, hgh + 1)

    f1 = fnprt(pg[1], ss; ylm = ylm, xlbl = "", ylbl = L"N_{\mathrm{parents}}", fdr = fdr)
    annotate!(15, 40, text(L"\Delta F = 0.005", 6))
    # f2 = fnprt(pg[2], ss; ylm = ylm, xlbl = "", ylbl = "", fdr = fdr)
    # annotate!(15, 40, text(L"\Delta F = 0.00625", 6))
    # f3 = fnprt(pg[3], ss; ylm = ylm, xlbl = "", ylbl = "", fdr = fdr)
    # annotate!(15, 40, text(L"\Delta F = 0.0075", 6))
    f4 = fnprt(pg[4], ss; ylm = ylm, ylbl = "", xlbl = "Generation", fdr = fdr)
    annotate!(15, 40, text(L"\Delta F = 0.01", 6))
    plot(f1, f4, layout = (1, 2), size = (600, 300))
end

function figtbv(; bs = "pg", nrs = "05", chr = "c01")
    pg = gatherdata("/home/xijiang/workspace/xyBnG/rst/paper-1/$bs/$nrs/$chr")
    fdr = -parse(Int, nrs)

    ss = keys(pg[1])
    low, hgh = lmtof(pg, ss, :mtbv)
    ylm = (low - 1, hgh + 1)

    f1 = fmtbv(pg[1], ss; fdr = fdr, ylm = ylm, xlbl = "", ylbl = L"\overline{\mathrm{TBV}}")
    #f2 = fmtbv(pg[2], ss; ylm = ylm, xlbl = "", ylbl = "")
    #f3 = fmtbv(pg[3], ss; ylm = ylm, xlbl = "", ylbl = "")
    f4 = fmtbv(pg[4], ss; fdr = fdr, ylm = ylm, xlbl = "Generation", ylbl = "")
    plot(f1, f4, layout = (1, 2), size = (600, 300))
end

function figvrg(; bs = "pg", nrs = "05", chr = "c01")
    pg = gatherdata("/home/xijiang/workspace/xyBnG/rst/paper-1/$bs/$nrs/$chr")
    fdr = -parse(Int, nrs)

    ss = keys(pg[1])
    low, hgh = lmtof(pg, ss, :vtbv, :genicv)
    ylm = (low - .1, hgh + .1)

    f1 = fvarg(pg[1], ss; fdr = fdr, ylm = ylm, xlbl = "")
    f4 = fvarg(pg[4], ss; fdr = fdr, ylm = ylm, xlbl = "Generation", ylbl = "")
    plot(f1, f4, layout = (1, 2), size = (600, 300))
end

function figbdr(; bs = "pg", nrs = "05", chr = "c01")
    pg = gatherdata("/home/xijiang/workspace/xyBnG/rst/paper-1/$bs/$nrs/$chr")
    fdr = -parse(Int, nrs)

    ss = keys(pg[1])
    lms = [lmtof(pg, ss, :floor)...; lmtof(pg, ss, :ceiling)...]

    f1 = fibdry(pg[1], ss, lms; fdr = fdr)
    #f2 = fbdry(pg[2], ss; ylm = ylm, xlbl = "", ylbl = "")
    #f3 = fbdry(pg[3], ss; ylm = ylm, xlbl = "", ylbl = "")
    f4 = fbdry(pg[4], ss, lms; fdr = fdr)
    plot(f1, f4, layout = (1, 2), size = (600, 300))
end

function figibd(; bs = "pg", nrs = "05", chr = "c01")
    pg = gatherdata("/home/xijiang/workspace/xyBnG/rst/paper-1/$bs/$nrs/$chr")
    fdr = -parse(Int, nrs)

    ss = keys(pg[1])
    low, hgh = lmtof(pg, ss, :fibd, :fped, :fhet, :fdrift, :fhet2, :fdrift2)
    f1 = fibrd(pg[1], ss; fdr = fdr, ylm = (low - .1, hgh + .1), xlbl = "", ylbl = L"\Delta F")
end
