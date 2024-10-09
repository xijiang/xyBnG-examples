"""
This is a test on `2024-10-03:02`. Both the random selection and OCS stages uses
`plnb`, which is `Plan(50, 50, 200)`.

## requires
- quick-test.jl
- xyBnG.jl v1.3.2

## commit
- 
"""
function trial()
    # 2024-10-03
    # includet("quick-test.jl")
    # commit: 
    qtest(rst = "quick-test/2024-10-04", nrpt = 20, nchp = 1100, nrng = 5)
end
