function itpos2coord(idxpos, tmitps)
    return [itp(idxpos) for itp in tmitps]
end

function z2itpos(z, itpz)
    zs = itpz.coefs
    zclamped = clamp(z, minimum(zs), maximum(zs))
    for i=2:length(zs)
        if (zs[i-1] <= zclamped <= zs[i]) || (zs[i] <= zclamped <= zs[i-1])
            return i-1 + (zclamped - zs[i-1]) / (zs[i] - zs[i-1])
        end
    end
end

function z2tmcoords(z, allitps)
    itpos = [z2itpos(z, tmitps[3]) for tmitps in allitps]
    return [itpos2coord(p, tmitps) for (p,tmitps) in zip(itpos, allitps)]
end

function sidechaincentroid(r::PDBResidue) # ignoring mass for now
    scatoms = filter(a -> a.atom ∉ ["N","CA","C","O"], r.atoms)
    return length(scatoms) == 0 ? nothing : sum([a.coordinates for a in scatoms]) / length(scatoms)
end

function scvector(r::PDBResidue)
    cacoord = r.atoms[findfirst(a -> a.atom === "CA", r.atoms)].coordinates
    sccoord = sidechaincentroid(r)
    return sccoord === nothing ? [0.0,0.0,0.0] : sccoord - cacoord
end

function res_inside_hull(sccoord, tmcoords) # tmcoords is a list of 2D points in the z-plane of sccoord 
    h = jarvismarch(tmcoords)
    return insidehull(sccoord[1:2], h)
end

function tm_res_is_inward(r::PDBResidue, itps)
    scvec = scvector(r)
    (isnothing(scvec) || sum(scvec)) == 0 && return false
    sccoord = alphacarbon_coordinates(r) .+ scvec
    return res_inside_hull(sccoord, z2tmcoords(sccoord[3], itps))
end

"""
    inward_tm_residues(seq, tmidxs)

Return an array of boolean[] indicating which residues (of those specified by `tmidxs`) are inward-facing.

`tmidxs` is a vector (typically of length 7, with each entry corresponding to a transmembrane region) of ranges of residue indices.
"""
function inward_tm_residues(seq, tmidxs)
    tmcoords = [alphacarbon_coordinates_matrix(seq[tm]) for tm in tmidxs]
    itps = [[interpolate(tm[i,:], BSpline((i === 3 ? Linear() : Cubic()))) for i=1:3] for tm in tmcoords]
    return [[tm_res_is_inward(r, itps) for r in seq[tm]] for tm in tmidxs]
end


function ecl_res_is_inward(r::PDBResidue, topcenter)
    scvec = scvector(r)
    (isnothing(scvec) || sum(scvec)) == 0 && return false
    accoord = alphacarbon_coordinates(r) 
    return dot(topcenter .- accoord, scvec) > 0
end


"""
    inward_elc_residues(seq, tmidxs)

Return an array of boolean[] indicating which residues (of those specified by `eclidxs`) are inward-facing (i.e. downward toward the opening of the binding pocket).

`eclidxs` is a vector (with each entry corresponding to an extracellular loop) of ranges of residue indices.
"""
function inward_ecl_residues(seq, eclidxs)
    tm_top_idxs = filter(x -> 0 < x < length(seq), vcat([[ecl[1] - 1, ecl[end] + 1] for ecl in eclidxs]...))    
    topcenter = mean(alphacarbon_coordinates_matrix(seq[tm_top_idxs]); dims=2)
    return [[ecl_res_is_inward(r, topcenter) for r in seq[ecl]] for ecl in eclidxs]
end
