## Minimum-distance alignment

"""
    tform = align(fixedpos::AbstractMatrix{Float64}, moving::Chain, sm::SequenceMapping)
    tform = align(fixed::Chain, moving::Chain, sm::SequenceMapping)

Return a rotated and shifted version of `moving` so that the centroids of
residues `moving[sm]` have least mean square error deviation from positions
`fixedpos` or those of residues `fixed`. `fixed` or `fixedpos` should include
just the residues or positions you want to align to, but `moving` should be an
entire chain.
"""
function align(fixedpos::AbstractMatrix{Float64}, moving::StructureLike, sm::SequenceMapping; seqname=nothing, warn::Bool=true)
    length(sm) == size(fixedpos, 2) || throw(DimensionMismatch("reference has $(size(fixedpos, 2)) positions, but `sm` has $(length(sm))"))
    movres = moving[sm]
    keep = map(!isnothing, movres)
    if !all(keep) && warn
        @warn "missing $(length(keep)-sum(keep)) anchors for sequence $seqname"
    end
    movpos = residue_centroid_matrix(movres[keep])
    fixedpos = fixedpos[:,keep]
    return Transformation(movpos, fixedpos, sm[keep], findall(keep))
end
align(fixed, moving::StructureLike, sm::SequenceMapping; kwargs...) = align(residue_centroid_matrix(fixed), moving, sm; kwargs...)

"""
    mapping = map_closest(mapto::StructureLike, mapfrom::StructureLike)

Return a vector `mapping[i] = (j, distij)`, matching the `i`th residue
in `mapto` to the `j`th residue in `mapfrom` and reporting the
distance between them. The mapping minimizes the sum of distances.
`mapto` and `mapfrom` must already be aligned for this to be meaningful.

If `mapfrom` is shorter than `mapto`, some `j`s will be 0, indicating
a skipped residue in `mapto`.
"""
function map_closest(mapto::StructureLike, mapfrom::StructureLike)
    refcenters = residue_centroid_matrix(mapto)
    seqcenters = residue_centroid_matrix(mapfrom)
    return map_closest(refcenters, seqcenters)
end

function map_closest(refcenters::AbstractMatrix, seqcenters::AbstractMatrix)
    D = pairwise(Euclidean(), refcenters, seqcenters; dims=2)
    assignment, _ = hungarian(D)
    fillval = convert(eltype(D), NaN)
    return collect(zip(assignment, [j == 0 ? fillval : D[i, j] for (i, j) in enumerate(assignment)]))
end

"""
    tform = align_closest(mapto::StructLike, mapfrom::StructLike; Dthresh=5)
    tform = align_closest(coordsto, coordsfrom; Dthresh=5)

Return the rigid transformation best aligning `mapfrom` to `mapto`. The transformation is computed from
residues matched by `map_closest`, using only residues closer than `Dthresh`.

Because the mapping is determined by distance, this can only "tweak" an already-close alignment.
"""
align_closest(mapto, mapfrom; kwargs...) = align_closest(residue_centroid_matrix(mapto), mapfrom; kwargs...)
align_closest(coordsto::AbstractMatrix, mapfrom; kwargs...) = align_closest(coordsto, residue_centroid_matrix(mapfrom); kwargs...)
align_closest(coordsto::AbstractMatrix, coordsfrom::AbstractMatrix; kwargs...) =
    align_closest(coordsto, coordsfrom, map_closest(coordsto, coordsfrom); kwargs...)

align_closest(mapto, mapfrom, mapping; kwargs...) = align_closest(residue_centroid_matrix(mapto), mapfrom, mapping; kwargs...)
align_closest(coordsto::AbstractMatrix, mapfrom, mapping; kwargs...) =
    align_closest(coordsto, residue_centroid_matrix(mapfrom), mapping; kwargs...)
function align_closest(coordsto::AbstractMatrix, coordsfrom::AbstractMatrix, mapping; Dthresh=5)
    toidx, fromidx = Int[], Int[]
    for (i, (j, d)) in pairs(mapping)
        j == 0 && continue
        d > Dthresh && continue
        push!(toidx, i)
        push!(fromidx, j)
    end
    return Transformation(coordsto[:,toidx], coordsfrom[:,fromidx], toidx, fromidx)
end

## Align to membrane

"""
    tform = align_to_membrane(chain::ChainLike, tms; extracellular=true)

Compute the rigid transformation `tform` needed to align `chain` to the
membrane, given the transmembrane segments `tms` as residue-indexes (e.g.,
`[37:61, 74:96, ...]`). `extracelluar` should be true if the N-terminus of the
protein is extracellular (`chain[first(tms[1])]` is at the extracellular face),
and false otherwise.

`applytransform!(chain, tform)` (or the `model` that includes `chain`) will
re-orient `chain` so that the center of the membrane is `z=0` and extracellular
is positive.

The algorithm finds the membrane normal `u` by maximizing the ratio

    sumᵢ (u ⋅ vᵢ)²
  -----------------
    sumᵢ (u ⋅ δᵢ)²

where `vᵢ` is a vector parallel to the `i`th TM helix, and `δᵢ` is a
within-leaflet atomic displacement.
"""
function align_to_membrane(chain::ChainLike, tms; extracellular::Bool=true)
    function leaflet_displacements(residues)
        coords = reduce(hcat, [coordarray(r) for r in residues])
        meanpos = mean(coords, dims=2)
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
    B = delta_i*delta_i' + delta_e*delta_e'
    λ, u = eigen(A, B)
    u = u[:, argmax(λ)]
    u = u / norm(u)
    # Find the projection by `u` of atom coordinates in each leaflet
    proj_e = reduce(hcat, [coordarray(r) for r in leaflet_e])' * u
    proj_i = reduce(hcat, [coordarray(r) for r in leaflet_i])' * u
    # Compute a rotation that aligns u to the z-axis
    t = [0, 0, mean(proj_e) > mean(proj_i) ? 1 : -1]   # extracellular is positive
    q = QuatRotation(1 + u'*t, cross(u, t)...)
    # Apply the rotation to the coordinates, and offset the z-coordinate to the mean of the leaflets
    ccoords = coordarray(chain)
    μ = vec(mean(ccoords, dims=2))
    return Transformation(μ, [0, 0, (median(proj_e) + median(proj_i))/2 - μ'*u], RotMatrix(q))
end

"""
    tform = align_to_axes(strct)

Compute the transformation needed to align the principle axes of inertia of
`strct` with the coordinate axes.
"""
function align_to_axes(strct)
    scoords = coordarray(strct)
    μ = mean(scoords, dims=2)
    scoords = scoords .- μ
    rotmtrx = GaussianMixtureAlignment.inertial_transforms(scoords)[1]
    return Transformation(μ, [0, 0, 0], rotmtrx)
end

function collect_leaflet_residues(chain, tms, extracellular::Bool)
    leaflet_e, leaflet_i = Residue[], Residue[]
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

## Needleman-Wunsch alignment

# This implements sequence alignment based on three-dimensional structure, where
# the cost of aligning two residues is the Euclidean distance between their
# centroids. The Needleman-Wunsch algorithm is used to find the optimal alignment
# of two sequences.

@enum NWParent NONE DIAG LEFT UP
Base.show(io::IO, p::NWParent) = print(io, p == NONE ? 'X' :
                                           p == DIAG ? '↖' :
                                           p == LEFT ? '←' :
                                           p == UP   ? '↑' : '?')
Base.show(io::IO, ::MIME"text/plain", p::NWParent) = show(io, p)

struct NWGapCosts{T}
    extend1::T
    extend2::T
    # Gap-opening costs
    open1::T
    open2::T

    function NWGapCosts{T}(extend1, extend2, open1, open2) where T
        all(>=(zero(T)), (extend1, extend2, open1, open2)) || throw(ArgumentError("gap costs must be non-negative"))
        return new{T}(extend1, extend2, open1, open2)
    end
end

"""
    gapcosts = NWGapCosts{T}(; extend1=0, extend2=0, open1=0, open2=0)

Create an affine cost for gaps in Needleman-Wunsch alignment. The cost of a gap of
length `k` is

    extend * k + open

All costs must be nonnegative.

`gapcosts(ϕ, idxs1, idxs2)` computes the contribution of gaps to the cost of
alignment `ϕ` between two sequences with `idxs1 = eachindex(seq1)` and `idxs2 =
eachindex(seq2)`. (The indices are needed to determine whether the alignment
starts or ends with a gap.)
"""
NWGapCosts{T}(; extend1=0, extend2=0, open1=0, open2=0) where T =
    NWGapCosts{T}(extend1, extend2, open1, open2)

Base.eltype(::Type{<:NWGapCosts{T}}) where T = T
Base.eltype(gapcosts::NWGapCosts) = eltype(typeof(gapcosts))

function (gapcosts::NWGapCosts{T})(ϕ, idxs1, idxs2) where T
    isempty(ϕ) && return saturating_add2(saturating_add2(gapcosts.ends1, gapcosts.ends2), length(idxs1)*gapcosts.extend1 + length(idxs2)*gapcosts.extend2)
    cost = zero(T)
    ϕprev = (0, 0)
    for ϕk in vcat(ϕ, last.((idxs1, idxs2)))
        dϕ = ϕk .- ϕprev
        cost = saturating_add2(cost, saturating_add2((dϕ[1] > 1) * gapcosts.open2, max(dϕ[1]-1, 0) * gapcosts.extend2))
        cost = saturating_add2(cost, saturating_add2((dϕ[2] > 1) * gapcosts.open1, max(dϕ[2]-1, 0) * gapcosts.extend1))
        ϕprev = ϕk
    end
    return cost
end

"""
    ϕ = align_nw(D, gapcosts::NWGapCosts)

Given a pairwise penalty matrix `D` (e.g., a pairwise distance matrix) and costs
for opening and extending gaps, find the optimal pairings

    ϕ = [(i1, j1), (i2, j2), ...]

that minimize

    sum(D[ϕk...] for ϕk in ϕ) + gapcosts(ϕ, axes(D)...)

subject to the constraint that `all(ϕ[k+1] .> ϕ[k])` for all `k`.
"""
function align_nw(D::AbstractMatrix, gapcosts::NWGapCosts)
    S, P = score_nw(D, gapcosts)
    S[end,end] == typemax(eltype(S)) && return nothing
    return traceback_nw(P)
end

"""
    ϕ = align_nw(seq1, seq2, gapcosts::NWGapCosts; mode=:distance_orientation)

Find the optimal `ϕ` matching `seq1[ϕ[k][1]]` to `seq2[ϕ[k][2]]` for all
`k`. `mode` controls the computation of pairwise matching penalties, and can be
either `:distance` or `:distance_orientation`, where the latter adds any
mismatch in sidechain orientation to the distance penalty.

`seq1` and `seq2` must be aligned to each other in 3D space before calling this
function. See [`align`](@ref).
"""
function align_nw(seq1::ChainLike, seq2::ChainLike, gapcosts::NWGapCosts;
                  mode::Symbol=:distance_orientation)
    supported_modes = (:distance, :distance_orientation)
    mode ∈ supported_modes || throw(ArgumentError("supported modes are $(supported_modes), got $mode"))

    refcenters = residue_centroid_matrix(seq1)
    srccenters = residue_centroid_matrix(seq2)
    D = pairwise(Euclidean(), refcenters, srccenters; dims=2)
    if mode == :distance_orientation
        refsc = reduce(hcat, [scvector(r) for r in seq1])
        srcsc = reduce(hcat, [scvector(r) for r in seq2])
        D = sqrt.(D.^2 + pairwise(SqEuclidean(), refsc, srcsc; dims=2))
    end
    return align_nw(D, gapcosts)
end

"""
    seqtms = align_ranges(seq1, seq2, seq2ranges::AbstractVector{<:AbstractUnitRange})

Transfer `refranges`, a list of reside index spans in `seq2`, to `seq1`. `seq1` and
`seq2` must be spatially aligned, and the assignment is made by minimizing
inter-chain distance subject to the constraint of preserving sequence order.
"""
function align_ranges(seq1::ChainLike, seq2::AbstractVector{<:AbstractResidue}, refranges::AbstractVector{<:AbstractUnitRange}; kwargs...)
    anchoridxs = sizehint!(Int[], length(refranges)*2)
    for r in refranges
        push!(anchoridxs, first(r), last(r))
    end
    issorted(anchoridxs) || throw(ArgumentError("`refranges` must be strictly increasing spans, got $refranges"))
    ϕ = align_nw(seq1, seq2[anchoridxs], NWGapCosts{Float64}(open1=Inf); kwargs...)
    @assert last.(ϕ) == eachindex(anchoridxs)
    return [ϕ[i][1]:ϕ[i+1][1] for i in 1:2:length(ϕ)]
end
align_ranges(seq1::ChainLike, seq2::Chain, refranges::AbstractVector{<:AbstractUnitRange}; kwargs...) =
    align_ranges(seq1, collectresidues(seq2), refranges; kwargs...)

function score_nw(D::AbstractMatrix, gapcosts::NWGapCosts)
    Base.require_one_based_indexing(D)
    m, n = size(D)
    T = typeof(zero(eltype(D)) + zero(eltype(gapcosts)))
    S = OffsetMatrix{T}(undef, 0:m, 0:n)            # S[i, j] is the optimal score of align(seq1[1:i], seq2[1:j])
    P = OffsetMatrix{NWParent}(undef, 0:m, 0:n)     # parent (predecessor of position i,j)
    fill!(S, -1)
    # Note: UP implies we're advancing sequence 1, and thus sequence 2 is in a gap
    #       LEFT is the converse
    # Initialize the scoring
    S[0,0] = zero(T); P[0,0] = NONE
    S[1,0] = saturating_add2(gapcosts.open2, gapcosts.extend2); P[1:m,0] .= UP
    for i = 2:m
        S[i,0] = saturating_add2(S[i-1,0], gapcosts.extend2)
    end
    S[0,1] = saturating_add2(gapcosts.open1, gapcosts.extend1); P[0,1:n] .= LEFT
    for j = 2:n
        S[0,j] = saturating_add2(S[0,j-1], gapcosts.extend1)
    end
    # Compute the score to each midpoint position
    for i = 1:m, j = 1:n
        left = saturating_add2(S[i, j-1], gapcosts.extend1)
        if P[i, j-1] != LEFT
            left = saturating_add2(left, gapcosts.open1)
        end
        up = saturating_add2(S[i-1,j], gapcosts.extend2)
        if P[i-1, j] != UP
            up = saturating_add2(up, gapcosts.open2)
        end
        diag = saturating_add2(S[i-1,j-1], D[i, j])
        minscore = min(diag, left, up)
        P[i, j] = minscore == diag ? DIAG :
                  minscore == left ? LEFT : UP
        S[i, j] = minscore
    end
    return S, P
end

function traceback_nw(P; start=last(CartesianIndices(P)))
    # Trace the path
    ϕ = Tuple{Int,Int}[]
    i, j = Tuple(start)
    istop, jstop = first.(axes(P))
    while i > istop && j > jstop
        p = P[i, j]
        if p == DIAG
            push!(ϕ, (i, j))
            i -= 1
            j -= 1
        elseif p == LEFT
            j -= 1
        elseif p == UP
            i -= 1
        else
            break
        end
    end
    return reverse!(ϕ)
end

# Saturating arithmetic to prevent overflow of Integer types
# https://discourse.julialang.org/t/the-performance-of-saturating-operations-or-adding-intrinsics/48575
# (thanks!)
saturating_add2(x::T, y::T) where {T <: Unsigned} = x + min(~x, y)

function saturating_add2(x::T, y::T) where {T <: Signed}
	x + ifelse(x < zero(x), max(y, typemin(x) - x), min(y, typemax(x) - x))
end

saturating_add2(x::T, y::T) where {T <: AbstractFloat} = x + y

saturating_add2(x::Real, y::Real) = saturating_add2(promote(x, y)...)
