# Plots for paper I

using DataFrames
using LaTeXStrings
using Measures
using Plots
using Serialization
using Statistics
"""
readRst()
In the directory which has a structure like this:
```
├── c01-1
│   ├── 1
│   └── 4
├── c01-2
│   ├── 1
│   └── 4
└── c29
    ├── 1
    └── 4
```
This function reads the summary files, merge the first 2, and returns a
    dictionary of DataFrames.
    """
function readRst()
    r1 = deserialize("c01-1/1/summary.ser")
    r2 = deserialize("c01-2/1/summary.ser")
    r2.repeat .+= 50
    append!(r1, r2)

    r2 = deserialize("c01-1/4/summary.ser")
    r3 = deserialize("c01-2/4/summary.ser")
    r3.repeat .+= 50
    append!(r2, r3)

    r3 = deserialize("c29/1/summary.ser")
    r4 = deserialize("c29/4/summary.ser")

    mpg, vpg = [], []
    for r in (r1, r2, r3, r4)
        mps = Dict{String,DataFrame}()
        vps = Dict{String,DataFrame}()
        nrpt = r.repeat[end]
        for grp in groupby(r, :scheme)
            mps[grp.scheme[1]] = combine(
                groupby(grp, :grt),
                Not(:repeat, :scheme) .=> mean .=> Not(:repeat, :scheme)
                )
            tdf = combine(
                groupby(grp, :grt),
                Not(:repeat, :scheme, :grt, :nid) .=> std .=> Not(:repeat, :scheme, :grt, :nid)
                )
            vps[grp.scheme[1]] = select(tdf, Not(:grt)) ./ sqrt(nrpt)
        end
        push!(mpg, mps)
        push!(vpg, vps)
    end

    mpg, vpg # mean and mean standard deviation
end

function fig_mtbv(mpg, vpg)
    xtrm, n = Float64[], length(mpg)
    for i in 1:n
        for s in keys(mpg[i])
            df = mpg[i][s]
            append!(xtrm, extrema(df.mtbv))
        end
    end
    ylm = begin
        t = extrema(xtrm)
        d = t[2] - t[1]
        (t[1] - 0.02d, t[2] + 0.02d)
    end
    fs = []
    for i in 1:n
        fig = plot(
            #xlabel = "Generation",
            #ylabel = L"\overline{\mathrm{MTBV}}",
            legend = false,
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        for s in keys(mpg[i])
            df = mpg[i][s]
            plot!(fig, df.grt .- 5, df.mtbv, label = uppercase(s[1:2]))
        end
        push!(fs, fig)
    end
    plot(fs..., layout = (1, 2), size = (800, 400))
    savefig("mtbv.pdf")
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