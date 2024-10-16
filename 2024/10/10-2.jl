"""
    sumvchr(; vdir = "xyBnG/rst/var-chr/")
Run ``Plots.scalefontsizes(0.3)`` once to scale the font sizes.
"""
function sumvchr(; vdir = "xyBnG/rst/var-chr/")
    rst = Any[]
    nchr = [1:3; 4:5:29]
    for i in nchr
        dir = joinpath(vdir, string(i))
        push!(rst, getdata(dir))
    end
    ss = keys(rst[1])
    Plots.scalefontsizes(0.3) # only run once
    nprt = Any[]
    for dic in rst
        push!(nprt, fnprt(dic, ss, ylm = (13, 42)))
    end
    plot(nprt..., layout = (3, 3))
    #=
    savefig("/home/xijiang/docs/reseach-notes/assets/2024/10/10/nprt.png")
    
    mtbv = Any[]
    for dic in rst
        push!(mtbv, fmtbv(dic, ss, ylm = (-0.4, 15)))
    end
    plot(mtbv..., layout = (3, 3))
    savefig("/home/xijiang/docs/reseach-notes/assets/2024/10/10/mtbv.png")
    bdry = Any[]
    for dic in rst
        push!(bdry, fbdry(dic, ss))
    end
    plot(bdry..., layout = (3, 3))
    =#
end
