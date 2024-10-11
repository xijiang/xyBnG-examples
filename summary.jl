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
function getdata(dir)
    rst = deserialize(joinpath(dir, "summary.ser"))
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
    fbdry(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Boundaries")
Plot the Boundaries of each scheme in `ss`. The solid lines are lower boundaries. 
The dashed lines are upper boundaries.
"""
function fbdry(mps, ss; fdr = -5, ylm = :best, xlbl = "Generation", ylbl = "Boundaries")
    fig = plot(dpi = 300, ylim = ylm, xlabel = xlbl, ylabel = ylbl, legend = false)
    color = 1
    for s in ss
        df = mps[s]
        plot!(fig, df.grt .+ fdr, df.floor, color = color)
        annotate!(df.grt[end] .+ fdr, df.floor[end], text(uppercase(s[1:2]), 10, :top))
        plot!(fig, df.grt .+ fdr, df.ceiling, color = color, linestyle = :dash)
        annotate!(df.grt[end] .+ fdr, df.ceiling[end], text(uppercase(s[1:2]), 10, :bottom))
        color += 1
    end
    fig
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
        annotate!(df.grt[end] .+ fdr, df.vtbv[end], text(uppercase(s[1:2]), 10, :top))
        plot!(fig, df.grt .+ fdr, df.genicv, linestyle = :dash, color = color)
        annotate!(df.grt[end] .+ fdr, df.genicv[end], text(uppercase(s[1:2]), 10, :bottom))
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
