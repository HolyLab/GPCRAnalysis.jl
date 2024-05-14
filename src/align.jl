"""
    newchain = align_to_membrane(chain::AbstractVector{PDBResidue}, tms; extracellular=true)

Align the `chain` to the membrane, given the transmembrane segments `tms` as
residue-indexes (e.g., `[37:61, 74:96, ...]`). `extracelluar` should be true if
the N-terminus of the protein is extracellular (`chain[first(tms[1])]` is at the
extracellular face), and false otherwise.

`newchain` is oriented so that the center of the membrane is `z=0` and
extracellular is positive.

The algorithm finds the membrane normal `u` by maximizing the ratio

    sumᵢ (u ⋅ vᵢ)²
  -----------------
    sumᵢ (u ⋅ δᵢ)²

where `vᵢ` is a vector parallel to the `i`th TM helix, and `δᵢ` is a within-leaflet
atomic displacement.
"""
function align_to_membrane(chain::AbstractVector{PDBResidue}, tms; extracellular::Bool=true)
    function leaflet_displacements(residues)
        coords = reduce(vcat, [coordinatesmatrix(r) for r in residues])
        meanpos = mean(coords, dims=1)
        return coords .- meanpos
    end

    # Collect the atoms at each leaflet face among the transmembrane segments
    leaflet_e, leaflet_i = collect_leaflet_residues(chain, tms, extracellular)
    # Collect the atom displacements in each leaflet
    delta_e = leaflet_displacements(leaflet_e)
    delta_i = leaflet_displacements(leaflet_i)
    # Also collect the vectors for the TM helices
    tm_vectors = [residue_centroid(chain[first(r)]) - residue_centroid(chain[last(r)]) for r in tms]
    # Find the membrane normal u, maximizing the ratio
    #         sumᵢ (u ⋅ vᵢ)²
    #       -----------------
    #         sumᵢ (u ⋅ δᵢ)²
    # where vᵢ is a vector of the TM helix and δᵢ is a within-leaflet displacement vector.
    # This is a generalized eigenvalue problem, which can be solved by finding the
    # eigenvector corresponding to the smallest eigenvalue of the matrix.
    μtm = mean(tm_vectors)
    A = sum(u -> (Δu = u - μtm; Δu * Δu'), tm_vectors)
    B = delta_i'*delta_i + delta_e'*delta_e
    λ, u = eigen(A, B)
    u = u[:, argmax(λ)]
    u = u / norm(u)
    # Find the projection by `u` of atom coordinates in each leaflet
    proj_e = reduce(vcat, [coordinatesmatrix(r) for r in leaflet_e]) * u
    proj_i = reduce(vcat, [coordinatesmatrix(r) for r in leaflet_i]) * u
    # Compute a rotation that aligns u to the z-axis
    t = [0, 0, mean(proj_e) > mean(proj_i) ? 1 : -1]   # extracellular is positive
    q = QuatRotation(1 + u'*t, cross(u, t)...)
    # Apply the rotation to the coordinates, and offset the z-coordinate to the mean of the leaflets
    coords = coordinatesmatrix(chain)
    μ = vec(mean(coords, dims=1))
    newcoords = (q * (coords' .- μ))' .+ [0, 0, (median(proj_e) + median(proj_i))/2 - μ'*u]'
    return change_coordinates(chain, newcoords)
end

function align_to_axes(strct)
    coords = coordinatesmatrix(strct)
    coords = coords .- mean(coords, dims=1)
    tform = GaussianMixtureAlignment.inertial_transforms(coords')[1]
    tform = AffineMap(Matrix(tform.linear),[tform.translation...])
    return change_coordinates(strct, tform(coords')')
end

function collect_leaflet_residues(chain, tms, extracellular)
    leaflet_e, leaflet_i = PDBResidue[], PDBResidue[]
    for r in tms
        start, stop = extrema(r)
        if extracellular
            push!(leaflet_e, chain[start])
            push!(leaflet_i, chain[stop])
        else
            push!(leaflet_i, chain[start])
            push!(leaflet_e, chain[stop])
        end
        extracellular = !extracellular
    end
    return leaflet_e, leaflet_i
end
