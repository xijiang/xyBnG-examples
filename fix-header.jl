using xyBnG
"""
    fix_header(file)
Fix the header of a XY file, if it is not correct. This happens when collecting
tskit or other base population data. The genotype matrix is stored as a BitArray,
but the header is not updated.
"""
function fix_header(file)
    hdr = xyBnG.XY.header(file)
    nlc, nhp = xyBnG.XY.dim(file)
    fz = filesize(file)
    if hdr.type ≠ 13
        sz = sizeof(xyBnG.xyTypes._type(hdr.type))
        if fz == sz * nlc * nhp + 24
            @info "This is a proper file, no fix needed"
            return
        end
        if abs(fz - nlc * nhp / 8 - 24) < 10
            @info "This is maybe of BitArray, header type will be changed ..."
            hdr.r, hdr.type = 1, 13
            xyBnG.XY.header!(file, hdr)
            @info "Done"
        end
    else
        if abs(fz - nlc * nhp / 8 - 24) < 10
            @info "This is a proper file, no fix needed"
        else
            @info "I don't know how to fix this file"
        end
    end
end
