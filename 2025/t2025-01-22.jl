"""
Results for paper I.

- quick test with empirical F0 for A-, G-, and I-OCSS
- 1 chromosome only
- ΔF = 1.0%, 0.5%

Uses paper-i.jl.
"""
function T20250122()
    pge([1.58534110], 0.005, 5, "rst/pge/c01/1")
    pge([1.58534110], 0.01,  5, "rst/pge/c01/4")
end
