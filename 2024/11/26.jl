using DataFrames
using LaTeXStrings
using Measures
using Plots
using Serialization
using Statistics

#=
Questions to answer:
	1. What are the differences between the forward and backward simulator
	2. How different \Delta\overline{\mathrm{TBV}} of different schemes are
	3. What about the genetic potential left?
	4. Which inbreeding indicators are to be used? How different they are.
	5. What are the recommendations.
	6. Different inbreeding restrictions, what are the differences between one chromosome and 29 chromosomes.
	7. Differences about different methods in translating variances into progresses.
Plotting:
	1. Plot the shades, or confident intervals.

## NOTE: include this file only.
=#

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

function lmtof(pg, ss, par...; fg = 1)
    ys = Float64[]
    for s in ss
        df = pg[s]
        for p in par
            append!(ys, df[!, p][fg:end])
        end
    end
    extrema(ys)
end

"""
Plot ceiling, ΔTBV, and floor in tight layout.
"""
function ctf(; file = "summary.jl")
    nrs = 5
    pg = getdata("."; file = file)
    ss = keys(pg)
    low, hgh = lmtof(pg, ss, :ceiling)
    ylm = (low, hgh) .* (.99, 1.01)
    h = Float64[]

    # Ceiling
    fclg = plot(
        dpi = 300,
        ylim = ylm,
        legend = false,
        xaxis = false,
        ytickfontsize = 6,
    )
    for s in ss
        df = pg[s]
        plot!(fclg, df.grt .- nrs, df.ceiling)
    end
    push!(h, hgh - low)
    # ΔTBV
    low, hgh = lmtof(pg, ss, :mtbv)
    ftbv = plot(
        dpi = 300,
        ylim = (low - 1, hgh + 1),
        legend = false,
        xaxis = false,
        ytickfontsize = 6,
    )
    for s in ss
        df = pg[s]
        plot!(ftbv, df.grt .- nrs, df.mtbv)
    end
    push!(h, hgh - low)
    # Floor
    low, hgh = lmtof(pg, ss, :floor)
    fflr = plot(dpi = 300, ylim = (low - 1, hgh + 1), ytickfontsize = 6, left_margin = 5mm)
    for s in ss
        df = pg[s]
        plot!(fflr, df.grt .- nrs, df.floor, label = uppercase(s[1:2]))
    end
    push!(h, hgh - low)
    fclg, ftbv, fflr, h
end

function mtctf(da, db, dc, dd)
    f1, f2, f3, hs = ctf(file = da)
    annotate!(f1, (-0.2, 0.9), text("Ceiling", 10, :right, rotation = 90))
    annotate!(f2, (-0.2, 0.5), text("ΔTBV", 10, :top, rotation = 90))
    annotate!(f3, (-0.2, 0.1), text("Floor", 10, :left, rotation = 90))
    fa = plot(f1, f2, f3, layout = grid(3, 1, heights = hs / sum(hs)), size = (400, 600))
    f1, f2, f3, hs = ctf(file = db)
    fb = plot(f1, f2, f3, layout = grid(3, 1, heights = hs / sum(hs)), size = (400, 600))
    f1, f2, f3, hs = ctf(file = dc)
    fc = plot(f1, f2, f3, layout = grid(3, 1, heights = hs / sum(hs)), size = (400, 600))
    f1, f2, f3, hs = ctf(file = dd)
    annotate!(f3, (0.8, 0.1), text("Generation", 10, :top))
    fd = plot(f1, f2, f3, layout = grid(3, 1, heights = hs / sum(hs)), size = (400, 600))
    plot(fa, fb, fc, fd, layout = (1, 4), size = (1200, 600))
end

function nprt(fdr, nrs, nch; k = 1)
    dir = joinpath(fdr, nrs, nch)
    pg = gatherdata(dir)
    ss = keys(pg[k])
    sft = parse(Int, nrs) + 2

    low, hgh = lmtof(pg, ss, :nsire, :ndam; fg = sft)
    ylm = (low - 1, hgh + 1)
    fig = plot(dpi = 300, ylim = ylm)
    clr = 1
    for s in ss
        df = pg[k][s]
        plot!(
            fig,
            df.grt[sft:end] .- sft .+ 1,
            df.nsire[sft:end],
            label = uppercase(s[1:2]),
            color = clr,
        )
        plot!(
            fig,
            df.grt[sft:end] .- sft .+ 1,
            df.ndam[sft:end],
            linestyle = :dot,
            color = clr,
            label = false,
        )
        clr += 1
    end
    fig
end

function mtprt(fdr, nrs)
    f1 = nprt(fdr, nrs, "c01")
    f2 = nprt(fdr, nrs, "c29")
    f3 = nprt(fdr, nrs, "c01", k = 4)
    f4 = nprt(fdr, nrs, "c29", k = 4)
    plot(f1, f2, f3, f4, layout = (1, 4), size = (1200, 300))
end

function vbv(fdr, nrs, nch; k = 1)
    dir = joinpath(fdr, nrs, nch)
    pg = gatherdata(dir)
    ss = keys(pg[k])
    sft = parse(Int, nrs) + 1

    low, hgh = lmtof(pg, ss, :vtbv, :genicv)
    ylm = (low - .1, hgh + .1)
    gtc = plot(dpi = 300, ylim = ylm)
    gnc = plot(dpi = 300, ylim = ylm)
    for s in ss
        df = pg[k][s]
        plot!(
            gtc,
            df.grt .- sft,
            df.vtbv,
            label = uppercase(s[1:2]),
        )
        plot!(
            gnc,
            df.grt .- sft,
            df.genicv,
            label = false,
        )
    end
    plot(gtc, gnc, layout = (2, 1), size = (400, 600))
end

function mtvbv(fdr, nrs)
    f1 = vbv(fdr, nrs, "c01")
    f2 = vbv(fdr, nrs, "c29")
    f3 = vbv(fdr, nrs, "c01", k = 4)
    f4 = vbv(fdr, nrs, "c29", k = 4)
    plot(f1, f2, f3, f4, layout = (1, 4), size = (1200, 600))
end

function ibd(fdr, nrs, nch; k =1)
    dir = joinpath(fdr, nrs, nch)
    pg = gatherdata(dir)
    ss = keys(pg[k])
    sft = parse(Int, nrs) + 1

    low, hgh = lmtof(pg, ss, :fibd, :fdrift2, :fhet2)
    ylm = (low, hgh) .* (.99, 1.01)
    fibd = plot(dpi = 300, ylim = ylm)
    fhet = plot(dpi = 300, ylim = ylm)
    fdft = plot(dpi = 300, ylim = ylm)

    for s in ss
        df = pg[k][s]
        plot!(
            fibd,
            df.grt .- sft,
            df.fibd,
            label = uppercase(s[1:2]),
        )
        plot!(
            fdft,
            df.grt .- sft,
            df.fdrift2,
            label = false,
        )
        plot!(
            fhet,
            df.grt .- sft,
            df.fhet2,
            label = false,
        )
    end
    plot(fibd, fhet, fdft, layout = (1, 3), size = (400, 600))
end

function mtibd(fdr, nrs)
    f1 = ibd(fdr, nrs, "c01")
    f2 = ibd(fdr, nrs, "c29")
    f3 = ibd(fdr, nrs, "c01", k = 4)
    f4 = ibd(fdr, nrs, "c29", k = 4)
    plot(f1, f2, f3, f4, layout = (1, 4), size = (1200, 600))
end