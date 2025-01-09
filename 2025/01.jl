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

"""
    line_clr(dat)
Assign a fixed color to each line in `dat`. The colors are for all the figures.
"""
function line_clr(dat)
    cs = Dict{String, Int}()
    for (i, k) in enumerate(keys(dat))
        cs[k] = i
    end
    cs
end

"""
    xtrm(dat, col)
Get the extrema of the column `col` in all the DataFrames in `dat`. Also add a
2% margin to the extrema.
"""
function xtrm(dat, col)
    x, n = Float64[], length(dat)
    for i in 1:n
        for s in keys(dat[i])
            df = dat[i][s]
            append!(x, extrema(df[!, col]))
        end
    end
    ylm = begin
        t = extrema(x)
        d = t[2] - t[1]
        (t[1] - 0.02d, t[2] + 0.02d)
    end
    ylm
end

"""
    ord_last(data, col)
Order the legends of the last plot in `data` by the last value of `col` in
reverse order.
"""
function ord_last(data, col)
    ks = collect(keys(data))
    dt = zeros(length(ks))
    for (i, k) in enumerate(ks)
        df = data[k]
        dt[i] = df[!, col][end]
    end
    ks = ks[sortperm(dt, rev = true)]
end

"""
    fig_mtbv(mpg, vpg, clr)
Plot the mean total breeding value of each generation of each scheme in `mpg`.
The shaded areas are the standard errors.
Each scheme is plotted in a fixed color defined in `clr`.
"""
function fig_mtbv(mpg, vpg, clr)
    ylm = xtrm(mpg, :mtbv)
    fs, n = [], length(mpg)
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :mtbv)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.mtbv,
                fillalpha = 0.2,
                ribbon = vpg[i][s].mtbv,
                label = uppercase(s[1:2]),
                legend = :bottomright,
                foreground_color_legend = nothing,
                background_color_legend = nothing,
                color = clr[s])
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -4.5, ylm[2] * .9,
        text(L"\overline{\mathrm{TBV}}", 10, :left, rotation = 90))
    annotate!(fs[2], 15, 0, text("Generation", 10, :bottom))
    plot(fs..., layout = (1, 2), size = (800, 400))
    savefig("mtbv.pdf")
    savefig("mtbv.png")
end

"""
    fig_nprt(mpg, vpg, clr)
Plot the number of parents of each generation of each scheme in `mpg`. The
shaded areas are the standard errors. Each scheme is plotted in a fixed color
defined in `clr`. The solid lines are the number of sires and the dashed lines
are the number of dams. Legends are in two columns. They read column by column.
"""
function fig_nprt(mpg, vpg, clr; p = 0)
    ylm = begin
        t = collect(xtrm(mpg, :nsire))
        append!(t, xtrm(mpg, :ndam))
        extrema(t)
    end
    fs, n, rg = [], length(mpg), 7:35
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :nsire)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt[rg] .- 5,
                df.nsire[rg],
                fillalpha = 0.2,
                ribbon = vpg[i][s].nsire,
                label = uppercase(s[1:2]),
                foreground_color_legend = nothing,
                background_color_legend = nothing,
                legend = abs(p - i) == 2 ? :right : :topright,
                color = clr[s])
            plot!(
                fig,
                df.grt[rg] .- 5,
                df.ndam[rg],
                linestyle = :dash,
                label = false,
                color = clr[s])
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        1.5, ylm[2] * .8,
        text(L"\overline{N_{\mathrm{parent}}}", 10, :left, rotation = 90))
    annotate!(fs[2], 25, 1, text("Generation", 10, :bottom))
    plot(fs..., layout = (1, 2), size = (800, 400))
    savefig("nprt.pdf")
    savefig("nprt.png")
end

"""
    fig_vg(mpg, vpg, clr)
Plot the genic and genetic variance of each generation of each scheme in `mpg`.
The genic variances were shaded with their mean standard deviation. Each scheme
is plotted in a fixed color defined in `clr`. The solid lines are the genic var
and the dashed lines are genetic variances. Legends are in the bottom left
corner.
"""
function fig_vg(mpg, vpg, clr)
    ylm = begin
        t = collect(xtrm(mpg, :vtbv))
        append!(t, xtrm(mpg, :genicv))
        extrema(t)
    end
    fs, n = [], length(mpg)
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :genicv)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.genicv,
                fillalpha = 0.1,
                ribbon = vpg[i][s].genicv,
                label = uppercase(s[1:2]),
                color = clr[s])
            plot!(
                fig,
                df.grt .- 5,
                df.vtbv,
                linestyle = :dash,
                label = false,
                foreground_color_legend = nothing,
                background_color_legend = nothing,
                legend = :topright,
                color = clr[s])
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -4.5,
        ylm[1] + 0.33 * (ylm[2] - ylm[1]),
        text("Genic/genetic variance", 10, :left, rotation = 90))
    h = (ylm[2] - ylm[1]) * 0.02 + ylm[1]
    annotate!(fs[2], 19, h, text("Generation", 10, :bottom))
    plot(fs..., layout = (1, 2), size = (800, 400))
    savefig("vg.pdf")
    savefig("vg.png")
end

"""
    fig_space(mpg, vpg, clr)
Plot the space limited by breeding ceiling and floor of each generation of each
scheme in `mpg`. The shaded areas are the mean standard deviations. Each scheme
is plotted in a fixed color defined in `clr`. The solid lines are the breeding
ceiling and the dashed lines are breeding floor. Legends are in the bottom left
corner.
"""
function fig_space(mpg, vpg, clr; p = -10)
    ylm = begin
        t = collect(xtrm(mpg, :ceiling))
        append!(t, xtrm(mpg, :floor))
        extrema(t)
    end
    fs, n = [], length(mpg)
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :ceiling)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.ceiling,
                fillalpha = 0.2,
                ribbon = vpg[i][s].ceiling,
                label = uppercase(s[1:2]),
                legend = :left,
                color = clr[s])
            end
            ks = ord_last(mpg[i], :floor)
            for s in ks
                df = mpg[i][s]
                plot!(
                    fig,
                    df.grt .- 5,
                    df.floor,
                    linestyle = :dash,
                    fillalpha = 0.2,
                    ribbon = vpg[i][s].floor,
                    label = uppercase(s[1:2]),
                    color = clr[s])
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        p, ylm[2] * .57, #-8 for 3:4
        text("Brd. val.", 10, :left, rotation = 90))
    annotate!(fs[1], 20, ylm[1] * .52, text("Floor", 10, :top))
    annotate!(fs[1], 20, ylm[2] * .52, text("Ceiling", 10, :bottom))
    annotate!(fs[2], 25, ylm[1], text("Generation", 10, :bottom))
    plot(fs..., layout = (1, 2), size = (800, 400))
    savefig("space.pdf")
    savefig("space.png")
end

"""
    fig_inbreeding(mpg, vpg, clr)
Inbreeding:
    - IBD: Identity by descent
    - Pedigree: Pedigree based inbreeding
    - Het.: Heterozygosity
    - Drift: Genetic drift
Plot the inbreeding of each generation of each scheme in `mpg`. The shaded
areas are the mean standard deviations. Each scheme is plotted in a fixed color
defined in `clr`. Legends are in the top left corner.
"""
function fig_inbreeding(mpg, vpg, clr)
    ylm = begin
        t = collect(xtrm(mpg, :fibd))
        append!(t, xtrm(mpg, :fped))
        append!(t, xtrm(mpg, :fhet3))
        append!(t, xtrm(mpg, :fdrift3))
        extrema(t)
    end
    fs, n = [], length(mpg)
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :fibd)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fibd,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fibd,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s])
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("IBD", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :fped)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fped,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fped,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s])
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Pedigree", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :fhet3)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fhet3,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fhet3,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s])
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Het.", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :fdrift3)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fdrift3,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fdrift3,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s])
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Drift", 10))
        push!(fs, fig)
    end
    h = (ylm[2] - ylm[1]) * 0.05 + ylm[1]
    annotate!(fs[8], 21, h, text("Generation", 8, :bottom))
    plot(fs..., layout = (2, 4), size = (800, 400))
    savefig("inbreeding.pdf")
    savefig("inbreeding.png")
end

"""
    fig_fxqtl(mpg, vpg, clr)
Plot the number of fixed QTL of each generation of each scheme in `mpg`. The
shaded areas are the mean standard deviations. Each scheme is plotted in a fixed
color defined in `clr`. The solid lines are the number of favorite QTL and the
dashed lines are the number of unfavorite QTL. Legends are in the right side.
"""
function fig_fxqtl(mpg, vpg, clr)
    ylm = begin
        t = collect(xtrm(mpg, :xfq))
        append!(t, xtrm(mpg, :xuq))
        extrema(t)
    end
    fs, n = [], length(mpg)
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :xfq)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.xfq,
                fillalpha = 0.2,
                ribbon = vpg[i][s].xfq,
                label = uppercase(s[1:2]),
                legendfontsize = 6,
                legend = :right,
                color = clr[s])
        end
        annotate!(fig, 15, ylm[2] * .2, text("Favorite", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :xuq)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.xuq,
                linestyle = :dash,
                fillalpha = 0.2,
                ribbon = vpg[i][s].xuq,
                label = uppercase(s[1:2]),
                legendfontsize = 6,
                legend = :right,
                color = clr[s])
        end
        annotate!(fig, 15, ylm[2] * .2, text("Un-f...e", 10))
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -3, ylm[2] * .65,
        text("Number of fixed QTL", 8, :left, rotation = 90))
    annotate!(fs[4], 20, 0, text("Generation", 8, :bottom))
    plot(fs..., layout = (1, 4), size = (800, 400))
    savefig("fixed-qtl.pdf")
    savefig("fixed-qtl.png")
end

"""
    fig_prp_fixed(mpg, vpg, clr; s = 1)
Plot the proportion of loci fixed of each generation of each scheme in `mpg`. The
shaded areas are the mean standard deviations. Each scheme is plotted in a fixed
color defined in `clr`. The solid lines are the proportion of QTL fixed, the
dashed lines are the proportion of reference fixed, and the dotted lines are the
proportion of chip fixed. Legends are in the bottom right corner.
"""
function fig_prp_fixed(mpg, vpg, clr; s = 1)
    tt = s == 1 ? (637, 3184) : (10000, 50000)
    ylm = begin
        t = collect(xtrm(mpg, :xqtl)) ./ tt[1]
        append!(t, xtrm(mpg, :xref) ./ tt[1])
        append!(t, xtrm(mpg, :xchip) ./ tt[2])
        extrema(t)
    end
    fs, n = [], length(mpg)
    
    for i in 1:n
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :xqtl)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.xqtl ./ tt[1],
                fillalpha = 0.2,
                ribbon = vpg[i][s].xqtl ./ tt[1],
                label = uppercase(s[1:2]),
                legend = :bottomright,
                legendfontsize = 6,
                color = clr[s])
        end
        annotate!(fig, 15, ylm[2] * .4, text("QTL", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :xref)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.xref ./ tt[1],
                fillalpha = 0.2,
                ribbon = vpg[i][s].xref ./ tt[1],
                label = uppercase(s[1:2]),
                legend = :bottomright,
                legendfontsize = 6,
                color = clr[s])
        end
        annotate!(fig, 15, ylm[2] * .4, text("Ref.", 10))
        push!(fs, fig)
        fig = plot(
            dpi = 300,
            ylim = ylm,
            size = (600, 400),
        )
        ks = ord_last(mpg[i], :xchip)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.xchip ./ tt[2],
                fillalpha = 0.2,
                ribbon = vpg[i][s].xchip ./ tt[2],
                label = uppercase(s[1:2]),
                legend = :bottomright,
                legendfontsize = 6,
                color = clr[s])
        end
        annotate!(fig, 15, ylm[2] * .4, text("Chip", 10))
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -3.5, ylm[2] * .45,
        text("Prop. loci fixed", 8, :left, rotation = 90))
    annotate!(fs[6], 15, 0, text("Generation", 8, :bottom))
    plot(fs..., layout = (2, 3), size = (800, 400))
    savefig("prp-fixed-$s.pdf")
    savefig("prp-fixed-$s.png")
end

"""
    fig_frq_evo()
Plot the evolution of allele frequency changes during the course of selection.
The MAF are categorized into 399 groups. That is the occurrence number of of
segregating allele 1. This is normalized by the total number of segregating
loci. Only the frequencies of generation 1 and 30 are plotted.
"""
function fig_frq_evo(clr)
    for k in keys(clr)
    end
end

"""
    all_fig()
Generate all the figures for the paper.
"""
function all_fig()
    mpg, vpg = readRst()
    clr = line_clr(mpg[1])

    begin # Figure of mean TBV
        fig_mtbv(mpg[1:2], vpg[1:2], clr)
        mv("mtbv.pdf", "mtbv-01.pdf", force = true)
        mv("mtbv.png", "mtbv-01.png", force = true)
        fig_mtbv(mpg[3:4], vpg[3:4], clr)
        mv("mtbv.pdf", "mtbv-29.pdf", force = true)
        mv("mtbv.png", "mtbv-29.png", force = true)
    end

    begin # Figure of number of parents
        fig_nprt(mpg[1:2], vpg[1:2], clr)
        mv("nprt.pdf", "nprt-01.pdf", force = true)
        mv("nprt.png", "nprt-01.png", force = true)
        fig_nprt(mpg[3:4], vpg[3:4], clr; p = 3)
        mv("nprt.pdf", "nprt-29.pdf", force = true)
        mv("nprt.png", "nprt-29.png", force = true)
    end

    begin # Figure of genic and genetic variance
        fig_vg(mpg[1:2], vpg[1:2], clr)
        mv("vg.pdf", "vg-01.pdf", force = true)
        mv("vg.png", "vg-01.png", force = true)
        fig_vg(mpg[3:4], vpg[3:4], clr)
        mv("vg.pdf", "vg-29.pdf", force = true)
        mv("vg.png", "vg-29.png", force = true)
    end

    begin # Figure of space between breeding ceiling and floor
        fig_space(mpg[1:2], vpg[1:2], clr)
        mv("space.pdf", "space-01.pdf", force = true)
        mv("space.png", "space-01.png", force = true)
        fig_space(mpg[3:4], vpg[3:4], clr; p = -8)
        mv("space.pdf", "space-29.pdf", force = true)
        mv("space.png", "space-29.png", force = true)
    end

    begin # Figure of inbreeding
        fig_inbreeding(mpg[1:2], vpg[1:2], clr)
        mv("inbreeding.pdf", "inbreeding-01.pdf", force = true)
        mv("inbreeding.png", "inbreeding-01.png", force = true)
        fig_inbreeding(mpg[3:4], vpg[3:4], clr)
        mv("inbreeding.pdf", "inbreeding-29.pdf", force = true)
        mv("inbreeding.png", "inbreeding-29.png", force = true)
    end

    begin # Figure of fixed QTL
        fig_fxqtl(mpg[1:2], vpg[1:2], clr)
        mv("fixed-qtl.pdf", "fixed-qtl-01.pdf", force = true)
        mv("fixed-qtl.png", "fixed-qtl-01.png", force = true)
        fig_fxqtl(mpg[3:4], vpg[3:4], clr)
        mv("fixed-qtl.pdf", "fixed-qtl-29.pdf", force = true)
        mv("fixed-qtl.png", "fixed-qtl-29.png", force = true)
    end

    begin # Figure of propotion loci fixed
        fig_prp_fixed(mpg[1:2], vpg[1:2], clr)
        mv("prp-fixed-1.pdf", "prp-fixed-01.pdf", force = true)
        mv("prp-fixed-1.png", "prp-fixed-01.png", force = true)
        fig_prp_fixed(mpg[3:4], vpg[3:4], clr; s = 2)
        mv("prp-fixed-2.pdf", "prp-fixed-29.pdf", force = true)
        mv("prp-fixed-2.png", "prp-fixed-29.png", force = true)
    end
end
