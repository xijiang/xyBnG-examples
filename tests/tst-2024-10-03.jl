"""
This is a test on `2024-10-03`

## requires
- quick-test.jl
- modified xyBnG.jl v1.3.2 with constraint = 2dF

## commit
"""
function devel()
    # 2024-10-03
    # includet("ggocs.jl")
    # commit: 
    qtest(rst = "quick-test/2df", nrpt = 20, nchp=1100)
end
