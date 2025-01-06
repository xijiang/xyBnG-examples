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
            vps[grp.scheme[1]] = select(tdf, Not(:grt)) ./ nrpt
        end
        push!(mpg, mps)
        push!(vpg, vps)
    end

    mpg, vpg # mean and mean standard deviation
end

function the_plots()
end
