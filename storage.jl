"""
    gb_storage(nrepeat, nscheme, ngrt, nlc, nid)
Calculate the storage requirement for `xyBnG` requirement if mid-results are
kept.
"""
function gb_storage(nrepeat, nscheme, ngrt, nlc, nid)
    # unique coding is 4 bytes
    nid * nlc * 4 * ngrt * nscheme * nrepeat / 1024^3
end
