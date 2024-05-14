## Minimum-distance alignment

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
    R = MIToS.PDB.kabsch((fixedpos .- refmean)', (movpos .- movmean)')
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

## Align to membrane

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

## One-sided Needleman-Wunsch alignment

@enum NWParent NONE DIAG LEFT
Base.show(io::IO, p::NWParent) = print(io, p == NONE ? 'X' :
                                           p == DIAG ? '↖' :
                                           p == LEFT ? '←' : '?')
Base.show(io::IO, ::MIME"text/plain", p::NWParent) = show(io, p)

"""
    ϕ = align_nw(D)

Given a penalty matrix `D` (e.g., a pairwise distance matrix), find
the optimal `ϕ` that minimizes the sum of `D[i, ϕ[i]]` for all `i`, subject
to the constraint that `ϕ[i+1] > ϕ[i]`. No gaps in `i` are allowed.
We require `size(D, 1) <= size(D, 2)`.
"""
function align_nw(D::AbstractMatrix)
    m, n = size(D)
    S, P = score_nw(D)
    # Trace the path
    ϕ = Vector{Int}(undef, m)
    i, j = m, n
    while i > 0
        P[i,j] == DIAG && (ϕ[i] = j; i -= 1; j -= 1)
        P[i,j] == LEFT && (j -= 1)
        P[i,j] == NONE && break
    end
    return ϕ
end

"""
    ϕ = align_nw(refseq, srcseq; mode=:distance)

Find the optimal `ϕ` matching `refseq[i]` to `srcseq[ϕ[i]]` for all `i`, subject
to the constraint that `ϕ[i+1] > ϕ[i]`.
"""
function align_nw(refseq::AbstractVector{PDBResidue}, srcseq::AbstractVector{PDBResidue}; mode::Symbol=Symbol("distance+scvector"))
    supported_modes = (:distance, Symbol("distance+scvector"))
    mode ∈ supported_modes || throw(ArgumentError("supported modes are $(supported_modes), got $mode"))

    refcenters = residue_centroid_matrix(refseq)
    srccenters = residue_centroid_matrix(srcseq)
    D = pairwise(Euclidean(), refcenters, srccenters; dims=2)
    if mode == Symbol("distance+scvector")
        refsc = reduce(hcat, [scvector(r) for r in refseq])
        srcsc = reduce(hcat, [scvector(r) for r in srcseq])
        D = sqrt.(D.^2 + pairwise(SqEuclidean(), refsc, srcsc; dims=2))
    end
    return align_nw(D)
end

"""
    seqtms = align_ranges(seq, ref, refranges::AbstractVector{<:AbstractUnitRange})

Transfer `refranges`, a list of reside index spans in `ref`, to `seq`. `seq` and
`ref` must be spatially aligned, and the assignment is made by minimizing
inter-chain distance subject to the constraint of preserving sequence order.
"""
function align_ranges(seq::AbstractVector{PDBResidue}, ref::AbstractVector{PDBResidue}, refranges::AbstractVector{<:AbstractUnitRange}; kwargs...)
    anchoridxs = sizehint!(Int[], length(refranges)*2)
    for r in refranges
        push!(anchoridxs, first(r), last(r))
    end
    issorted(anchoridxs) || throw(ArgumentError("`refranges` must be strictly increasing spans, got $refranges"))
    ϕ = align_nw(ref[anchoridxs], seq; kwargs...)
    return [ϕ[i]:ϕ[i+1] for i in 1:2:length(ϕ)]
end

function score_nw(D::AbstractMatrix{T}) where T
    Base.require_one_based_indexing(D)
    m, n = size(D)
    m <= n || throw(ArgumentError("First dimension cannot be longer than the second"))
    S = OffsetMatrix{T}(undef, 0:m, 0:n)
    P = OffsetMatrix{NWParent}(undef, 0:m, 0:n)
    # Initialize the scoring. We need to use all rows of D while always incrementing the columns.
    # Give maximum badness to paths that don't do this.
    S[:,0] .= typemax(T); S[0,:] .= zero(T)
    P[:,0] .= NONE; P[0,:] .= NONE
    for i = 1:m
        for j = 1:i-1
            S[i,j], P[i,j] = typemax(T), NONE
        end
        for j = n-m+i+1:n
            S[i,j], P[i,j] = typemax(T), NONE
        end
    end
    # Compute the score of each possible path
    for i = 1:m, j = i:n-m+i
        d = D[i,j]
        sd, sl = S[i-1,j-1], S[i, j-1]
        if sd == typemax(T)
            S[i,j], P[i,j] = sl, LEFT
        else
            diag = sd + d
            left = sl
            S[i,j], P[i,j] = diag < left ? (diag, DIAG) : (left, LEFT)
        end
    end
    return S, P
end
