# Plots for paper I

using DataFrames
using LaTeXStrings
using Mmap
using Measures
using Plots
using Serialization
using Statistics

function mrgRst(clr, dA, dB, dC)
    mkpath.(["$dC/1", "$dC/4"]; mode = 0o700)
    ks = collect(keys(clr))
    append!(ks, ["founder", "rand"])
    for r in (1, 4)
        rst, nrpt, w = begin
            t1 = deserialize("$dA/$r/summary.ser")
            n1 = t1.repeat[end]
            t2 = deserialize("$dB/$r/summary.ser")
            n2 = t2.repeat[end]
            t2.repeat .+= n1
            append!(t1, t2)
            t1, Dict(dA => (0, n1), dB => (n1, n2)), ndigits(n1 + n2)
        end
        serialize("$dC/$r/summary.ser", rst)
        
        for (dir, (p, n)) in nrpt
            v = ndigits(n)
            for i in 1:n
                t1 = lpad(i, v, '0')
                t2 = lpad(i + p, w, '0')
                for k in ks
                    mv("$dir/$r/$t1-$k.xy", "$dC/$r/$t2-$k.xy", force = true)
                    mv("$dir/$r/$t1-$k.ped", "$dC/$r/$t2-$k.ped", force = true)
                end
                mv("$dir/$r/$t1-snp.xy", "$dC/$r/$t2-snp.xy", force = true)
                mv("$dir/$r/$t1-founder.lmp", "$dC/$r/$t2-founder.lmp", force = true)
            end
        end
    end
end

"""
    fltRst(dir, clr)
Filter the results in `dir` by the schemes in `clr`. The genotype files in
`*.xy` are deleted. Make sure that the frequencies are summarized.
"""
function fltRst(dir, clr)
end

function readSup(dir)
    r1 = deserialize(joinpath(dir, "1/summary.ser"))
    r2 = deserialize(joinpath(dir, "4/summary.ser"))

    mpg, vpg = [], []
    for r in (r1, r2)
        mps = Dict{String,DataFrame}()
        vps = Dict{String,DataFrame}()
        nrpt = r.repeat[end]
        for grp in groupby(r, :scheme)
            mps[grp.scheme[1]] = combine(
                groupby(grp, :grt),
                Not(:repeat, :scheme) .=> mean .=> Not(:repeat, :scheme),
            )
            tdf = combine(
                groupby(grp, :grt),
                Not(:repeat, :scheme, :grt, :nid) .=>
                    std .=> Not(:repeat, :scheme, :grt, :nid),
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
    cs = Dict{String,Int}()
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
    for i = 1:n
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
    fig_mtbv(mpg, vpg, clr; c = 1)
Plot the mean total breeding value of each generation of each scheme in `mpg`.
The shaded areas are the standard errors.
Each scheme is plotted in a fixed color defined in `clr`.
"""
function fig_mtbv(mpg, vpg, clr; c = 1)
    ylm = xtrm(mpg, :mtbv)
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -4.5,
        ylm[2] * 0.9,
        text(L"\overline{\mathrm{TBV}}", 10, :left, rotation = 90),
    )
    p = ylm[2] * 0.9
    annotate!(fs[1], 2, p, text("A: " * ch * L"F_{0.5\%}", 6))
    annotate!(fs[2], 15, 0, text("Generation", 10, :bottom))
    annotate!(fs[2], 2, p, text("B: " * ch * L"F_{1\%}", 6))
    plot(fs..., layout = (1, 2), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("mtbv-$chr.pdf")
    savefig("mtbv-$chr.png")
end

"""
    fig_nprt(mpg, vpg, clr; p = 0, c = 1)
Plot the number of parents of each generation of each scheme in `mpg`. The
shaded areas are the standard errors. Each scheme is plotted in a fixed color
defined in `clr`. The solid lines are the number of sires and the dashed lines
are the number of dams. Legends are in two columns. They read column by column.
"""
function fig_nprt(mpg, vpg, clr; p = 0, c = 1)
    ylm = begin
        t = collect(xtrm(mpg, :nsire))
        append!(t, xtrm(mpg, :ndam))
        extrema(t)
    end
    fs, n, rg = [], length(mpg), 7:35
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
            plot!(
                fig,
                df.grt[rg] .- 5,
                df.ndam[rg],
                linestyle = :dash,
                label = false,
                color = clr[s],
            )
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        1.5,
        ylm[2] * 0.8,
        text(L"\overline{N_{\mathrm{parent}}}", 10, :left, rotation = 90),
    )
    p = ylm[2] * 0.9
    annotate!(fs[1], 5, p, text("A: " * ch * L"F_{0.5\%}", 6))
    annotate!(fs[2], 25, 1, text("Generation", 10, :bottom))
    annotate!(fs[2], 5, p, text("B: " * ch * L"F_{1\%}", 6))
    plot(fs..., layout = (1, 2), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("nprt-$chr.pdf")
    savefig("nprt-$chr.png")
end

"""
    fig_vg(mpg, vpg, clr)
Plot the genic and genetic variance of each generation of each scheme in `mpg`.
The genic variances were shaded with their mean standard deviation. Each scheme
is plotted in a fixed color defined in `clr`. The solid lines are the genic var
and the dashed lines are genetic variances. Legends are in the bottom left
corner.
"""
function fig_vg(mpg, vpg, clr; c = 1)
    ylm = begin
        t = collect(xtrm(mpg, :vtbv))
        append!(t, xtrm(mpg, :genicv))
        extrema(t)
    end
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
            plot!(
                fig,
                df.grt .- 5,
                df.vtbv,
                linestyle = :dash,
                label = false,
                foreground_color_legend = nothing,
                background_color_legend = nothing,
                legend = :topright,
                color = clr[s],
            )
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -4.5,
        ylm[1] + 0.33 * (ylm[2] - ylm[1]),
        text("Genic/genetic variance", 10, :left, rotation = 90),
    )
    p = ylm[2] * 0.95
    annotate!(fs[1], 10, p, text("A: " * ch * L"F_{0.5\%}", 6))
    h = (ylm[2] - ylm[1]) * 0.02 + ylm[1]
    annotate!(fs[2], 19, h, text("Generation", 10, :bottom))
    annotate!(fs[2], 10, p, text("B: " * ch * L"F_{1\%}", 6))
    plot(fs..., layout = (1, 2), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("vg-$chr.pdf")
    savefig("vg-$chr.png")
end

"""
    fig_space(mpg, vpg, clr; p = -10, c = 1)
Plot the space limited by breeding ceiling and floor of each generation of each
scheme in `mpg`. The shaded areas are the mean standard deviations. Each scheme
is plotted in a fixed color defined in `clr`. The solid lines are the breeding
ceiling and the dashed lines are breeding floor. Legends are in the bottom left
corner.
"""
function fig_space(mpg, vpg, clr; p = -10, c = 1)
    ylm = begin
        t = collect(xtrm(mpg, :ceiling))
        append!(t, xtrm(mpg, :floor))
        extrema(t)
    end
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
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
                color = clr[s],
            )
        end
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        p,
        ylm[2] * 0.57, #-8 for 3:4
        text("Brd. val.", 10, :left, rotation = 90),
    )
    pp = ylm[2] * 0.9
    annotate!(fs[1], 5, pp, text("A: " * ch * L"F_{0.5\%}", 6))
    annotate!(fs[1], 20, ylm[1] * 0.52, text("Floor", 10, :top))
    annotate!(fs[1], 20, ylm[2] * 0.52, text("Ceiling", 10, :bottom))
    annotate!(fs[2], 25, ylm[1], text("Generation", 10, :bottom))
    annotate!(fs[2], 5, pp, text("B: " * ch * L"F_{1\%}", 6))
    plot(fs..., layout = (1, 2), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("space-$chr.pdf")
    savefig("space-$chr.png")
end

"""
    fig_inbreeding(mpg, vpg, clr; c = 1)
Inbreeding:
    - IBD: Identity by descent
    - Pedigree: Pedigree based inbreeding
    - Het.: Heterozygosity
    - Drift: Genetic drift
Plot the inbreeding of each generation of each scheme in `mpg`. The shaded
areas are the mean standard deviations. Each scheme is plotted in a fixed color
defined in `clr`. Legends are in the top left corner.
"""
function fig_inbreeding(mpg, vpg, clr; c = 1)
    ylm = begin
        t = collect(xtrm(mpg, :aF))
        append!(t, xtrm(mpg, :iF))
        append!(t, xtrm(mpg, :iF2))
        append!(t, xtrm(mpg, :fhet2))
        append!(t, xtrm(mpg, :fdrift2))
        append!(t, xtrm(mpg, :fhet3))
        append!(t, xtrm(mpg, :fdrift3))
        extrema(t)
    end
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
        ks = ord_last(mpg[i], :iF)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.iF,
                fillalpha = 0.2,
                ribbon = vpg[i][s].iF,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("IBD-chip", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
        ks = ord_last(mpg[i], :iF2)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.iF2,
                fillalpha = 0.2,
                ribbon = vpg[i][s].iF2,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("IBD-ref", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
        ks = ord_last(mpg[i], :aF)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.aF,
                fillalpha = 0.2,
                ribbon = vpg[i][s].aF,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Pedigree", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
        ks = ord_last(mpg[i], :fhet2)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fhet2,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fhet2,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Het. ref", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Het. chip", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
        ks = ord_last(mpg[i], :fdrift2)
        for s in ks
            df = mpg[i][s]
            plot!(
                fig,
                df.grt .- 5,
                df.fdrift2,
                fillalpha = 0.2,
                ribbon = vpg[i][s].fdrift2,
                label = uppercase(s[1:2]),
                legend = :topleft,
                legendfontsize = 6,
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Drift ref", 10))
        push!(fs, fig)

        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        i == 2 && annotate!(fig, 15, ylm[2], text("Drift chip", 10))
        push!(fs, fig)
    end
    h = (ylm[2] - ylm[1]) * 0.05 + ylm[1]
    annotate!(fs[8], 21, h, text("Generation", 8, :bottom))
    p = ylm[2] * 0.85
    CS = 'A':'Z'
    for i in 1:7
        annotate!(fs[i], 15, p, text(CS[i] * ": " * ch * L"F_{0.5\%}", 6))
    end
    for i in 8:14
        annotate!(fs[i], 15, p, text(CS[i] * ": " * ch * L"F_{1\%}", 6))
    end
    plot(fs..., layout = (2, 7), size = (1100, 550))
    chr = c == 1 ? "01" : "29"
    savefig("inbreeding-$chr.pdf")
    savefig("inbreeding-$chr.png")
end

"""
    fig_fxqtl(mpg, vpg, clr; c = 1)
Plot the number of fixed QTL of each generation of each scheme in `mpg`. The
shaded areas are the mean standard deviations. Each scheme is plotted in a fixed
color defined in `clr`. The solid lines are the number of favorite QTL and the
dashed lines are the number of unfavorite QTL. Legends are in the right side.
"""
function fig_fxqtl(mpg, vpg, clr; c = 1)
    ylm = begin
        t = collect(xtrm(mpg, :xfq))
        append!(t, xtrm(mpg, :xuq))
        extrema(t)
    end
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"
    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        annotate!(fig, 15, ylm[2] * 0.2, text("Favorite", 10))
        push!(fs, fig)
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        annotate!(fig, 15, ylm[2] * 0.2, text("Unfavorite", 10))
        push!(fs, fig)
    end
    annotate!(
        fs[1],
        -3,
        ylm[2] * 0.65,
        text("Number of fixed QTL", 8, :left, rotation = 90),
    )
    p = ylm[2] * 0.9
    annotate!(fs[1], 10, p, text("A: " * ch * L"F_{0.5\%}", 6))
    annotate!(fs[2], 10, p, text("B: " * ch * L"F_{0.5\%}", 6))
    annotate!(fs[4], 20, 0, text("Generation", 8, :bottom))
    annotate!(fs[3], 10, p, text("C: " * ch * L"F_{1\%}", 6))
    annotate!(fs[4], 10, p, text("D: " * ch * L"F_{1\%}", 6))
    plot(fs..., layout = (1, 4), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("fixed-qtl-$chr.pdf")
    savefig("fixed-qtl-$chr.png")
end

"""
    fig_prp_fixed(mpg, vpg, clr; s = 1)
Plot the proportion of loci fixed of each generation of each scheme in `mpg`. The
shaded areas are the mean standard deviations. Each scheme is plotted in a fixed
color defined in `clr`. The solid lines are the proportion of QTL fixed, the
dashed lines are the proportion of reference fixed, and the dotted lines are the
proportion of chip fixed. Legends are in the bottom right corner.
"""
function fig_prp_fixed(mpg, vpg, clr; c = 1)
    tt = c == 1 ? (637, 3184) : (10000, 50000)
    ylm = begin
        t = collect(xtrm(mpg, :xqtl)) ./ tt[1]
        append!(t, xtrm(mpg, :xref) ./ tt[1])
        append!(t, xtrm(mpg, :xchip) ./ tt[2])
        extrema(t)
    end
    fs, n = [], length(mpg)
    ch = c == 1 ? L"C_{1}" : L"C_{29}"

    for i = 1:n
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        annotate!(fig, 15, ylm[2] * 0.4, text("QTL", 10))
        push!(fs, fig)
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        annotate!(fig, 15, ylm[2] * 0.4, text("Ref.", 10))
        push!(fs, fig)
        fig = plot(dpi = 300, ylim = ylm, size = (600, 400))
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
                color = clr[s],
            )
        end
        annotate!(fig, 15, ylm[2] * 0.4, text("Chip", 10))
        push!(fs, fig)
    end
    annotate!(fs[1], -3.5, ylm[2] * 0.45, text("Prop. loci fixed", 8, :left, rotation = 90))
    p = ylm[2] * 0.9
    CS = 'A':'F'
    for i in 1:3
        annotate!(fs[i], 5, p, text(CS[i] * ": " * ch * L"F_{0.5\%}", 6))
    end
    annotate!(fs[6], 15, 0, text("Generation", 8, :bottom))
    for i in 4:6
        annotate!(fs[i], 5, p, text(CS[i] * ": " * ch * L"F_{1\%}", 6))
    end
    plot(fs..., layout = (2, 3), size = (800, 400))
    chr = c == 1 ? "01" : "29"
    savefig("prp-fixed-$chr.pdf")
    savefig("prp-fixed-$chr.png")
end

"""
    sum_frq_evo(clr, ds, file)
Summarize the evolution of allele frequency changes during the course of
selection. The MAF are categorized into 399 groups. That is the occurrence
number of of segregating allele 1. This is normalized by the total number of
segregating loci. Only the frequencies of generation 1 and 30 are plotted.
"""
function sum_frq_evo(clr, dir, file)
    frq, ti, nhp = zeros(399, 2), zeros(Int, 3), 14400
    open(file, "w") do bin
        for r in (1, 4)
            rsm = deserialize("$dir/$r/summary.ser")
            nrpt = rsm.repeat[end]
            for k in keys(clr)
                print(' ', k)
                for i = 1:nrpt
                    tag = lpad(i, ndigits(nrpt), '0')
                    nlc = begin
                        read!("$dir/$r/$tag-$k.xy", ti)
                        ti[2]
                    end
                    xy = mmap("$dir/$r/$tag-$k.xy", Matrix{Int16}, (nlc, nhp), 24)
                    for s in sum(isodd.(xy[:, 2001:2400]), dims = 2)
                        (s == 0 || s == 400) && continue
                        frq[s, 1] += 1
                    end
                    for s in sum(isodd.(xy[:, 14001:14400]), dims = 2)
                        (s == 0 || s == 400) && continue
                        frq[s, 2] += 1
                    end
                end
                write(bin, frq)
                frq .= 0
            end
        end
    end
end

function moveavg(x, n)
    y = zeros(length(x) - n + 1)
    for i = 1:length(y)
        y[i] = mean(x[i:i+n-1])
    end
    y
end

function fig_frq_evo(clr, bin)
    frq = zeros(399, length(clr) * 4)
    cg1 = length(clr) + 1
    read!(bin, frq)
    fs = []
    ts = sum(frq[:, 1]) # Total segregating loci
    vs = vec(sum(frq[:, 2:2:12], dims = 1))
    ks = collect(keys(clr))[sortperm(vs, rev = true)]
    fig = plot(dpi = 300, size = (600, 400), legend = :top)
    ch = occursin("01", bin) ? L"C_{1}" : L"C_{29}"
    for k in ks
        scatter!(
            fig,
            1:399,
            frq[:, 2clr[k]] / ts,
            label = uppercase(k)[1:2],
            markerstrokewidth = 0,
            ms = 0.6,
            color = clr[k],
        )
        ma = moveavg(frq[:, 2clr[k]] / ts, 5)
        plot!(fig, 3:397, ma, color = clr[k], label = false)
    end
    scatter!(
        fig,
        1:399,
        frq[:, 1] / ts,
        label = L"G_1",
        markerstrokewidth = 0,
        ms = 0.6,
        color = cg1,
    )
    ma = moveavg(frq[1:399, 1] / ts, 5)
    plot!(fig, 3:397, ma, label = false, color = cg1)
    p = maximum(frq[:, 1] / ts) * 0.9
    annotate!(fig, 50, p, text("A: " * ch * L"F_{0.5\%}", 6))
    push!(fs, fig)
    ts = sum(frq[:, 13])
    vs = vec(sum(frq[:, 14:2:24], dims = 1))
    ks = collect(keys(clr))[sortperm(vs, rev = true)]
    fig = plot(dpi = 300, size = (600, 400), legend = :top)
    for k in ks
        scatter!(
            fig,
            1:399,
            frq[:, 2clr[k]+12] / ts,
            label = uppercase(k)[1:2],
            markerstrokewidth = 0,
            ms = 0.6,
            color = clr[k],
        )
        ma = moveavg(frq[:, 2clr[k]+12] / ts, 5)
        plot!(fig, 3:397, ma, color = clr[k], label = false)
    end
    scatter!(
        fig,
        1:399,
        frq[:, 13] / ts,
        label = L"G_{1}",
        markerstrokewidth = 0,
        ms = 0.6,
        color = cg1,
    )
    ma = moveavg(frq[:, 13] / ts, 5)
    plot!(fig, 3:397, ma, label = false, color = cg1)
    p = maximum(frq[:, 13] / ts) * 0.9
    annotate!(fig, 50, p, text("B: " * ch * L"F_{1\%}", 6))
    push!(fs, fig)
    plot(fs..., layout = (1, 2), size = (800, 400))
    chr = occursin("01", bin) ? "01" : "29"
    savefig("frq-evo-$chr.pdf")
    savefig("frq-evo-$chr.png")
end

function sup_fig(dir)
    mpg, vpg = readSup(dir)
    clr = line_clr(mpg[1]) # Assign a fixed color to each scheme
    chr = dir[end-1:end]
    c = chr == "01" ? 1 : 2

    fig_mtbv(mpg, vpg, clr; c = c)
    
    begin # Figure of number of parents
        p = chr == "01" ? 0 : 3
        fig_nprt(mpg, vpg, clr; p = p, c = c)
    end

    fig_vg(mpg, vpg, clr; c = c)

    begin # Figure of space between breeding ceiling and floor
        p = chr == "01" ? -10 : -8
        fig_space(mpg, vpg, clr; p = p, c = c)
    end

    fig_inbreeding(mpg, vpg, clr; c = c)

    fig_fxqtl(mpg, vpg, clr) # Figure of fixed QTL
    
    fig_prp_fixed(mpg, vpg, clr; c = c) # Figure of propotion loci fixed

    begin # Evolution of allele frequency
        isfile("ch-$chr.bin") || sum_frq_evo(clr, dir, "ch-$chr.bin")
        fig_frq_evo(clr, "ch-$chr.bin")
    end
end
