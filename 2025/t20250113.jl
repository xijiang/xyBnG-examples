"""
Results for paper I.

- quick test with empirical F0 for A-, G-, and I-OCSS
- 1 chromosome only
- ΔF = 1.0%, 0.5%

Uses paper-i.jl.
"""
function T2024112301(; nrng = 5, nrpt = 50)
    paper_1_pg([1.58534110], 0.005, nrng, "rst/pg/c01/1"; nrpt = nrpt)
    paper_1_pg([1.58534110], 0.01,  nrng, "rst/pg/c01/4"; nrpt = nrpt)
end
