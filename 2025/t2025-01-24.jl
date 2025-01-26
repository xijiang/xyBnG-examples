"""
Results for paper I.

- 1 vs 29 chromosomes
- pg (ulimit -s 1000000) vs tskit
- 5 vs 15 generations
- ΔF = 1.0%, 0.75%, 0.625%, 0.5%
- remove ri-, add ig-, not BLUPs    

Uses paper-i.jl.
"""
function Tamd0124()
    paper_1_pg([1.58534110], 0.005, 5, "rst/pgemp/c01/1"; nrpt = 100)
    paper_1_pg([1.58534110], 0.01,  5, "rst/pgemp/c01/4"; nrpt = 100)
end
