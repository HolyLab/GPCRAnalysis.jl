function itpos2coord(idxpos, tmitps)
    return [itp(idxpos) for itp in tmitps]
end

function z2itpos(z, itpz::Interpolations.BSplineInterpolation{<:Any, <:Any, <:Any, <:BSpline{<:Linear}, <:Any})
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

function sidechaincentroid(r::AbstractResidue) # ignoring mass for now
    scatoms = collectatoms(r, !backboneselector)
    return length(scatoms) == 0 ? nothing : sum(coords(a) for a in scatoms) / length(scatoms)
end

function scvector(r::AbstractResidue)
    cacoord = alphacarbon_coordinates(r)
    sccoord = sidechaincentroid(r)
    return sccoord === nothing ? zero(cacoord) : sccoord - cacoord
end

function res_inside_hull(sccoord, tmcoords) # tmcoords is a list of 2D points in the z-plane of sccoord
    h = jarvismarch(tmcoords)
    return insidehull(sccoord[1:2], h)
end

function tm_res_is_inward(r::AbstractResidue, itps)
    scvec = scvector(r)
    all(iszero, scvec) && return false
    sccoord = alphacarbon_coordinates(r) .+ scvec
    return res_inside_hull(sccoord, z2tmcoords(sccoord[3], itps))
end

"""
    inward_tm_residues(seq, tmidxs)

Return an array of boolean[] indicating which residues (of those specified by `tmidxs`) are inward-facing.

`tmidxs` is a vector (typically of length 7, with each entry corresponding to a transmembrane region) of ranges of residue indices.
"""
function inward_tm_residues(seq::AbstractVector{<:AbstractResidue}, tmidxs)
    tmcoords = [alphacarbon_coordinates_matrix(seq[tm]) for tm in tmidxs]
    itps = [[interpolate(tm[i,:], BSpline((i === 3 ? Linear() : Cubic()))) for i=1:3] for tm in tmcoords]
    return [[tm_res_is_inward(r, itps) for r in seq[tm]] for tm in tmidxs]
end
inward_tm_residues(seq::Chain, tmidxs) = inward_tm_residues(collectresidues(seq), tmidxs)

function ecl_res_is_inward(r::AbstractResidue, topcenter)
    scvec = scvector(r)
    all(iszero, scvec) && return false
    accoord = alphacarbon_coordinates(r)
    return dot(topcenter .- accoord, scvec) > 0
end


"""
    inward_ecl_residues(seq, eclidxs)

Return an array of boolean[] indicating which residues (of those specified by `eclidxs`) are inward-facing (i.e. downward toward the opening of the binding pocket).

`eclidxs` is a vector (with each entry corresponding to an extracellular loop) of ranges of residue indices.
"""
function inward_ecl_residues(seq::AbstractVector{<:AbstractResidue}, eclidxs; includes_amino_terminus::Bool=false)
    tm_top_idxs = filter(x -> 0 < x < length(seq), vcat([isempty(ecl) ? ecl : [ecl[begin] - 1, ecl[end] + 1] for ecl in eclidxs[begin+includes_amino_terminus:end]]...))
    topcenter = mean(alphacarbon_coordinates_matrix(seq[tm_top_idxs]); dims=2)
    return [[ecl_res_is_inward(r, topcenter) for r in seq[ecl]] for ecl in eclidxs]
end
inward_ecl_residues(seq, eclidxs; includes_amino_terminus::Bool=false) =
    inward_ecl_residues(collectresidues(seq), eclidxs; includes_amino_terminus)
