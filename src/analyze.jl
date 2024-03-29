"""
    X = project_sequences(msa::AbstractMultipleSequenceAlignment; fracvar::Real = 0.9)

Perform a classical multidimensional scaling analysis to project the sequences in `msa` to a space
in which pairwise distances approximately reproduce `100 - percentsimilarity(seq1, seq2)`.
The dimensionality is chosen to reconstruction `fracvar` of the variance.
"""
function project_sequences(msa::AbstractMultipleSequenceAlignment; fracvar::Real = 0.9)
    sim = percentsimilarity(msa)
    D = 100 .- Matrix(sim)
    f = fit(MDS, D; distances=true)
    # Capture sufficient variance
    cl = cumsum(f.λ)
    cl ./= cl[end]
    nd = findfirst(>=(fracvar), cl)
    X = predict(f)
    return X[1:nd, :]
end

const reduced_code = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

function _entropy(v)
    count = countitems(v)
    e = 0.0
    for (_, n) in count
        p = n / length(v)
        e += - p * log(p)
    end
    return e
end

"""
    columnwise_entropy(msa::AbstractMultipleSequenceAlignment, aacode = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

Compute the entropy of each column in an MSA. Low entropy indicates high conservation.

Unmatched entries (`'-'` residues) contribute to the entropy calculation as if they were an ordinary residue.
"""
function columnwise_entropy(msa::AbstractMultipleSequenceAlignment, aacode=reduced_code)
    resnum = map(r -> aacode[r], getresidues(msa))
    return map(_entropy, eachcol(resnum))
end

MSA.Residue(r::PDBResidue) = three2residue(r.id.name)

"""
    residue_centroid(r::PDBResidue)

Compute the mean position of all atoms in `r` excluding hydrogen.

Residue centers may yield a more reliable measure of "comparable residues" than the α-carbons (see `CAmatrix`) because
they incorporate the orientation of the side chain relative to the overall fold.

See also [`residue_centroid_matrix`](@ref).
"""
residue_centroid(r::PDBResidue) = mean(atom -> atom.coordinates, Iterators.filter(atom -> atom.element != "H", r.atoms))

"""
    residue_centroid_matrix(seq)

Return a matrix of all residue centroids as columns. See also [`residue_centroid`](@ref).
"""
residue_centroid_matrix(seq) = reduce(hcat, map(residue_centroid, seq))

"""
    align(fixedpos::AbstractMatrix{Float64}, moving::AbstractVector{PDBResidue}, sm::SequenceMapping)
    align(fixed::AbstractVector{PDBResidue}, moving::AbstractVector{PDBResidue}, sm::SequenceMapping)

Return a rotated and shifted version of `moving` so that the α-carbons of residues `moving[sm]` have least
mean square error deviation from positions `fixedpos` or those of residues `fixed`.
"""
function align(fixedpos::AbstractMatrix{Float64}, moving::AbstractVector{PDBResidue}, sm::SequenceMapping; seqname=nothing)
    length(sm) == size(fixedpos, 2) || throw(DimensionMismatch("reference has $(size(fixedpos, 2)) positions, but `sm` has $(length(sm))"))
    movres = moving[sm]
    keep = map(!isnothing, movres)
    if !all(keep)
        @warn "missing $(length(keep)-sum(keep)) anchors for sequence $seqname"
    end
    movpos = residue_centroid_matrix(movres[keep])
    fixedpos = fixedpos[:,keep]
    refmean = mean(fixedpos; dims=2)
    movmean = mean(movpos; dims=2)
    R = kabsch((fixedpos .- refmean)', (movpos .- movmean)')
    return change_coordinates(moving, (coordinatesmatrix(moving) .- movmean') * R .+ refmean')
end
align(fixed::AbstractVector{PDBResidue}, moving::AbstractVector{PDBResidue}, sm::SequenceMapping; kwargs...) =
    align(residue_centroid_matrix(fixed), moving, sm; kwargs...)

function mapclosest(mapto::AbstractVector{PDBResidue}, mapfrom::AbstractVector{PDBResidue})
    refcenters = residue_centroid_matrix(mapto)
    seqcenters = residue_centroid_matrix(mapfrom)
    D = pairwise(Euclidean(), refcenters, seqcenters; dims=2)
    assignment, _ = hungarian(D)
    return collect(zip(assignment, [j == 0 ? convert(eltype(D), NaN) : D[i, j] for (i, j) in enumerate(assignment)]))
    # d, m = findmin(D; dims=1)
    # return [mi[1] => di for (mi, di) in zip(vec(m), vec(d))]
end

"""
    chargelocations(pdb::AbstractVector{PDBResidue}; include_his::Bool=false)

Return a list of potential charge locations in the protein structure. Each is a tuple `(position, residueindex, AAname)`.
The positions are those of `N` (in positively-charged residues like Arg & Lys) or `O` (in negatively-charged residues like
Asp and Glu). N- and C-termini are not included in the list. While each residue will carry a net total charge of ±1,
the location of each potential charge will be listed (1 for Lys, 2 each for Arg, Asp, and Glu).

By default, histidine is not considered charged, but you can include it by setting `include_His=true`.
"""
function chargelocations(pdb::AbstractVector{PDBResidue}; include_His::Bool=false)
    argpat(atom) = startswith(atom.atom, "NH")  # arginine
    lyspat(atom) = atom.atom == "NZ"            # lysine
    asppat(atom) = startswith(atom.atom, "OD")  # aspartate
    glupat(atom) = startswith(atom.atom, "OE")  # glutamate
    hispat(atom) = startswith(atom.atom, r"N[DE]")  # histidine

    locs = Tuple{Coordinates,Int,String}[]
    for (i, res) in pairs(pdb)
        name = res.id.name
        if name == "ARG"
            idx = findfirst(argpat, res.atoms)
            push!(locs, (res.atoms[idx].coordinates, i, name))
            idx = findnext(argpat, res.atoms, idx+1)
            push!(locs, (res.atoms[idx].coordinates, i, name))
        elseif name == "LYS"
            idx = findfirst(lyspat, res.atoms)
            push!(locs, (res.atoms[idx].coordinates, i, name))
        elseif name == "ASP"
            idx = findfirst(asppat, res.atoms)
            push!(locs, (res.atoms[idx].coordinates, i, name))
            idx = findnext(asppat, res.atoms, idx+1)
            push!(locs, (res.atoms[idx].coordinates, i, name))
        elseif name == "GLU"
            idx = findfirst(glupat, res.atoms)
            push!(locs, (res.atoms[idx].coordinates, i, name))
            idx = findnext(glupat, res.atoms, idx+1)
            push!(locs, (res.atoms[idx].coordinates, i, name))
        elseif include_His && name == "HIS"
            idx = findfirst(hispat, res.atoms)
            push!(locs, (res.atoms[idx].coordinates, i, name))
            idx = findnext(hispat, res.atoms, idx+1)
            push!(locs, (res.atoms[idx].coordinates, i, name))
        end
    end
    return locs
end

positive_locations(chargelocs::AbstractVector{Tuple{Coordinates,Int,String}}) = [loc[1] for loc in chargelocs if loc[3] ∈ ("ARG", "LYS", "HIS")]
negative_locations(chargelocs::AbstractVector{Tuple{Coordinates,Int,String}}) = [loc[1] for loc in chargelocs if loc[3] ∈ ("ASP", "GLU")]
