using DataFrames
using Dates
using LaTeXStrings
using Plots
using Serialization
using Statistics

"""
    getdata(dir)
Read the summary file in the directory `dir` and return a dictionary of dataframes.
"""
function getdata(dir; file = "summary.ser")
    rst = deserialize(joinpath(dir, file))
    schemes = unique(rst.scheme)
    dic = Dict{String,DataFrame}()
    mps = Dict{String,DataFrame}()
    for scheme in schemes
        dic[scheme] = select(filter(row -> row.scheme == scheme, rst), Not(:scheme))
        mps[scheme] =
            combine(groupby(dic[scheme], :grt), Not(:repeat) .=> mean .=> Not(:repeat))
    end
    mps # mean of parameters in each scheme
end

function lmtof(data, ss, par...; fg = 1)
    ys = Float64[]
    for mps in data
        for s in ss
            df = mps[s]
            for p in par
                append!(ys, df[!, p][fg:end])
            end
        end
    end
    extrema(ys)
end

"""
    fnprt(mps, ss; fdr = -5, xlbl = "Generation", ylbl = "N_parents", ylm = :best)
Plot the number of parents in each generation of each scheme in `ss`. The solid
lines are number of sires. The dashed lines are number of dams.
"""
function fnprt(
    mps,
    ss;
    fdr = -5,
    xlbl = "Generation",
    ylbl = L"N_{\mathrm{parents}}",
    ylm = :best,
)
    fig = plot(xlabel = xlbl, ylabel = ylbl, legend = false, dpi = 300, ylim = ylm)
    color = 1
    for s in ss
        df = mps[s]
        rg = -fdr+2:length(df.grt)
        plot!(fig, df.grt[rg] .+ fdr, df.nsire[rg], linestyle = :solid, color = color)
        annotate!(df.grt[end] .+ fdr, df.nsire[end], text(uppercase(s[1:2]), 5, :top))
        plot!(fig, df.grt[rg] .+ fdr, df.ndam[rg], linestyle = :dot, color = color)
        color += 1
    end
    fig
end

"""
    fmtbv(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "mTBV")
Plot the mean TBV of each scheme in `ss`.
"""
function fmtbv(
    mps,
    ss;
    fdr = -5,
    ylm = :best,
    xlbl = "Generation",
    ylbl = L"\overline{\mathrm{TBV}}",
)
    fig = plot(dpi = 300, ylim = ylm, legend = false, xlabel = xlbl, ylabel = ylbl)
    for s in ss
        df = mps[s]
        plot!(fig, df.grt .+ fdr, df.mtbv)
        annotate!(df.grt[end] .+ fdr, df.mtbv[end], text(uppercase(s[1:2]), 4, :top))
    end
    fig
end

"""
    fbdry(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation")
Plot the Boundaries of each scheme in `ss`. The solid lines are lower boundaries. 
The dashed lines are upper boundaries.
"""
function fbdry(mps, ss, lms; fdr = -5)
    alm = (lms[3] - 5, lms[4] + 5) # ceilings
    blm = (lms[1] - 5, lms[2] + 5) # floors
    fc = plot(dpi = 300, ylim = alm, legend = false, size = (150, 300), xaxis = false)
    ff = plot(dpi = 300, ylim = blm, legend = false, size = (150, 300))
    color = 1
    for s in ss
        df = mps[s]
        plot!(
            ff,
            df.grt .+ fdr,
            df.floor,
            color = color,
        )
        annotate!(ff, df.grt[end] .+ fdr, df.floor[end], text(uppercase(s[1:2]), 5, :top))
        plot!(
            fc,
            df.grt .+ fdr,
            df.ceiling,
            color = color,
            linestyle = :dash,
        )
        annotate!(fc, df.grt[end] .+ fdr, df.ceiling[end], text(uppercase(s[1:2]), 5, :bottom))
        color += 1
    end
    plot(fc, ff, layout = (2, 1))
end

"""
    fvarg(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Var(TBV)")
Plot the variance of TBV of each scheme in `ss`. The solid lines are genetic
variances and the dashed lines are genic variances.
"""
function fvarg(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Var(TBV)")
    fig = plot(dpi = 300, ylim = ylm, legend = false, xlabel = xlbl, ylabel = ylbl)
    color = 1
    for s in ss
        df = mps[s]
        plot!(fig, df.grt .+ fdr, df.vtbv, color = color)
        annotate!(df.grt[end] .+ fdr, df.vtbv[end], text(uppercase(s[1:2]), 5, :top))
        plot!(fig, df.grt .+ fdr, df.genicv, linestyle = :dash, color = color)
        annotate!(df.grt[end] .+ fdr, df.genicv[end], text(uppercase(s[1:2]), 5, :bottom))
        color += 1
    end
    fig
end

"""
    fibrd(mps, ss, ps; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Inbreeding")
Plot the inbreeding of each scheme in `ss`. Line styles for `ps` in order of
[solid, dash, dot, dashdot].
"""
function fibrd(mps, ss, ps; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Inbreeding")
    fig = plot(dpi = 300, ylim = ylm, legend = false, xlabel = xlbl, ylabel = ylbl)
    color = 1
    ls = [:solid, :dash, :dot, :dashdot]
    for s in ss
        df, l = mps[s], 1
        for p in ps
            plot!(fig, df.grt .+ fdr, df[!, p], color = color, linestyle = ls[l])
            l += 1
        end
        annotate!(df.grt[end] .+ fdr, df[!, ps[1]][end], text(uppercase(s[1:2]), 10, :top))
        color += 1
    end
    fig
end

"""
    cmpibd(mps, s; fdr = -5, ylm = :best)
Compare inbreeding indicators of one scenario `s` in `mps`.
"""
function cmpibd(mps, s; fdr = -5, ylm = :best)
    fig = plot(dpi = 300, ylim = ylm, xlabel = "Generation", ylabel = "Inbreeding")
    plot!(mps[s].grt .+ fdr, mps[s].fibd, label = "IBD")
    annotate!(mps[s].grt[end] .+ fdr, mps[s].fibd[end], text("IBD", 5, :left))
    plot!(mps[s].grt .+ fdr, mps[s].fped, label = "Ped")
    annotate!(mps[s].grt[end] .+ fdr, mps[s].fped[end], text("Ped", 5, :left))
    # plot!(mps[s].grt .+ fdr, mps[s].fhet, label = "Het")
    # annotate!(mps[s].grt[end] .+ fdr, mps[s].fhet[end], text("Het", 5, :left))
    # plot!(mps[s].grt .+ fdr, mps[s].fdrift, label = "Drift")
    # annotate!(mps[s].grt[end] .+ fdr, mps[s].fdrift[end], text("Drift", 5, :left))
    plot!(mps[s].grt .+ fdr, mps[s].fhet2, label = "Het2")
    annotate!(mps[s].grt[end] .+ fdr, mps[s].fhet2[end], text("Het2", 5, :left))
    plot!(mps[s].grt .+ fdr, mps[s].fdrift2, label = "Drift2")
    annotate!(mps[s].grt[end] .+ fdr, mps[s].fdrift2[end], text("Drift2", 5, :left))
end
