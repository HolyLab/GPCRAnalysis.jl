"""
    X = project_sequences(msa; fracvar::Real = 0.9)

Perform a classical multidimensional scaling analysis to project the sequences in `msa` to a space
in which pairwise distances approximately reproduce `100 - percentsimilarity(seq1, seq2)`.
The dimensionality is chosen to reconstruction `fracvar` of the variance.
"""
function project_sequences(msa; fracvar::Real = 0.9)
    sim = percent_similarity(msa)
    D = 100 .- Matrix(sim)
    f = fit(MDS, D; distances=true)
    # Capture sufficient variance
    cl = cumsum(f.λ)
    cl ./= cl[end]
    nd = findfirst(>=(fracvar), cl)
    X = predict(f)
    return X[1:nd, :]
end

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
    columnwise_entropy(f, msa)

Compute the entropy of each column in an MSA, after applying `f` to each residue. Low entropy indicates high conservation.
"""
function columnwise_entropy(f, msa)
    resnum = map(f, residuematrix(msa))
    return map(_entropy, eachcol(resnum))
end

"""
    residue_centroid(r::AbstractResidue)

Compute the mean position of all atoms in `r` excluding hydrogen.

Residue centers may yield a more reliable measure of "comparable residues" than the α-carbons (see `CAmatrix`) because
they incorporate the orientation of the side chain relative to the overall fold.

See also [`residue_centroid_matrix`](@ref).
"""
residue_centroid(r::AbstractResidue) = vec(mean(coordarray(r, heavyatomselector); dims=2))

"""
    residue_centroid_matrix(seq)

Return a matrix of all residue centroids as columns. See also [`residue_centroid`](@ref).
"""
residue_centroid_matrix(seq) = reduce(hcat, map(residue_centroid, seq))

"""
    alphacarbon_coordinates(res::AbstractResidue)

Return the coordinates of the α-carbon in `res`.
"""
alphacarbon_coordinates(r) = coords(firstmatch(calphaselector, r))

"""
    alphacarbon_coordinates_matrix(seq)

Return a matrix of αC coordinates as columns across all residues. See also [`alphacarbon_coordinates`](@ref).
"""
alphacarbon_coordinates_matrix(seq) = reduce(hcat, map(alphacarbon_coordinates, seq))

"""
    chargelocations(chain::ChainLike; include_his::Bool=false)

Return a list of potential charge locations in the protein structure. Each is a tuple `(position, residueindex, AAname)`.
The positions are those of `N` (in positively-charged residues like Arg & Lys) or `O` (in negatively-charged residues like
Asp and Glu). N- and C-termini are not included in the list. While each residue will carry a net total charge of ±1,
the location of each potential charge will be listed (1 for Lys, 2 each for Arg, Asp, and Glu).

By default, histidine is not considered charged, but you can include it by setting `include_His=true`.
"""
function chargelocations(chain::ChainLike; include_His::Bool=false)
    argpat(atom) = startswith(atomname(atom), "NH")  # arginine
    lyspat(atom) = atomname(atom) == "NZ"            # lysine
    asppat(atom) = startswith(atomname(atom), "OD")  # aspartate
    glupat(atom) = startswith(atomname(atom), "OE")  # glutamate
    hispat(atom) = startswith(atomname(atom), r"N[DE]")  # histidine

    locs = Tuple{typeof(coords(first(first(chain)))),Int,String}[]
    for (i, res) in enumerate(chain)
        name = resname(res)
        if name == "ARG"
            for a in res
                argpat(a) && push!(locs, (coords(a), i, name))
            end
        elseif name == "LYS"
            for a in res
                lyspat(a) && push!(locs, (coords(a), i, name))
            end
        elseif name == "ASP"
            for a in res
                asppat(a) && push!(locs, (coords(a), i, name))
            end
        elseif name == "GLU"
            for a in res
                glupat(a) && push!(locs, (coords(a), i, name))
            end
        elseif include_His && name == "HIS"
            for a in res
                hispat(a) && push!(locs, (coords(a), i, name))
            end
        end
    end
    return locs
end

positive_locations(chargelocs::AbstractVector{<:Tuple{AbstractVector,Int,String}}) = [loc[1] for loc in chargelocs if loc[3] ∈ ("ARG", "LYS", "HIS")]
negative_locations(chargelocs::AbstractVector{<:Tuple{AbstractVector,Int,String}}) = [loc[1] for loc in chargelocs if loc[3] ∈ ("ASP", "GLU")]
