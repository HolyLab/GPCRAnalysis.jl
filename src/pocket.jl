function itpos2coord(idxpos, tmitps)
    return SVector(tmitps[1](idxpos), tmitps[2](idxpos))
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

function res_inside_hull(sccoord, tmcoords) # tmcoords is a list of SVector{2} points
    h = jarvismarch(tmcoords)
    return insidehull(SVector(sccoord[1], sccoord[2]), h)
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

function default_eclidxs(tmidxs)
    length(tmidxs) == 7 || error("Expected 7 TM helices in tmidxs, got $(length(tmidxs))")
    return (
        1:tmidxs[1][1]-1,
        tmidxs[2][end]+1:tmidxs[3][1]-1,
        tmidxs[4][end]+1:tmidxs[5][1]-1,
        tmidxs[6][end]+1:tmidxs[7][1]-1,
    )
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

function _pocket_z_range(seq::AbstractVector{<:AbstractResidue}, tmidxs)
    z_max = maximum(alphacarbon_coordinates(r)[3] for tm in tmidxs for r in seq[tm])
    return 0.0 .. z_max
end

"""
    cavity = detect_pocket_cavity(seq, tmidxs; kwargs...) -> IsotropicGMM{3,Float64}

Identify the binding pocket cavity of `seq` using a grid-based solvent-exclusion approach.

`seq` must be aligned to the membrane (see [`align_to_membrane`](@ref)) so that z is the
membrane normal and z > 0 is extracellular. Probe positions are placed on a regular grid
within the TM helix bundle and accepted if they do not clash with any protein atom.

Cavity detection has two stages:
1. **Geometric filter**: the probe must lie inside the 2D convex hull of the TM helix
   centerlines at that z-level (computed from `tmidxs` via B-spline interpolation).
2. **Clash filter**: the probe center must be at least `probe_radius + vdW_radius` away
   from every protein atom.

All 7 TM helices should normally be passed in `tmidxs` so that the convex hull correctly
encloses the full helix bundle.

# Keyword arguments
- `probe_radius=1.4`: solvent probe radius (Å); 1.4 approximates a water molecule
- `grid_spacing=0.25`: grid resolution (Å)
- `zi`: z-interval (Å) to search; defaults to the extracellular half of the TM bundle
  (`0 .. z_max` of the TM alpha carbons)
- `ρmax=12.0`: xy radius of the search cylinder (Å)
- `σ=1.4`: Gaussian σ assigned to each probe in the returned `IsotropicGMM`

The returned `IsotropicGMM` can be visualised with `gmmdisplay` (from
GaussianMixtureAlignment, requires a Makie backend) and used with `rocs_align` to
score steric complementarity of a ligand to the cavity.

See also: [`pocket_pharmacophore`](@ref), [`complementary_pharmacophore`](@ref).
"""
function detect_pocket_cavity(
    seq::AbstractVector{<:AbstractResidue},
    tmidxs;
    probe_radius::Real = 1.4,
    grid_spacing::Real = 0.25,
    zi::Union{AbstractInterval, Nothing} = nothing,
    ρmax::Real = 12.0,
    σ::Real = 1.4,
)
    zi === nothing && (zi = _pocket_z_range(seq, tmidxs))
    z_lo, z_hi = leftendpoint(zi), rightendpoint(zi)

    # Build TM centerline interpolants for the inside-bundle convex-hull test
    tmcoords_mats = [alphacarbon_coordinates_matrix(seq[tm]) for tm in tmidxs]
    itps = [[interpolate(tm[i,:], BSpline(i === 3 ? Linear() : Cubic())) for i=1:3]
            for tm in tmcoords_mats]

    # Collect atom coords and vdW radii for atoms near the search region
    buffer = probe_radius + 2.2   # ≥ max heavy-atom vdW radius
    atom_coords = Vector{SVector{3,Float64}}()
    atom_radii  = Vector{Float64}()
    for r in seq
        rname = resname(r)
        for a in r
            c = SVector{3,Float64}(coords(a))
            (hypot(c[1], c[2]) < ρmax + buffer &&
             z_lo - buffer ≤ c[3] ≤ z_hi + buffer) || continue
            aname = atomname(a) == "OXT" ? "O" : atomname(a)
            push!(atom_coords, c)
            push!(atom_radii, get(vanderwaalsradius, (rname, aname), 1.8))
        end
    end
    min_dist_sq = [(probe_radius + ar)^2 for ar in atom_radii]

    # Grid search: z is the outer loop so the convex hull is reused across the xy-plane
    gaussians = IsotropicGaussian{3,Float64}[]
    for z in range(z_lo, z_hi; step=Float64(grid_spacing))
        hull = jarvismarch(z2tmcoords(z, itps))
        for x in range(-ρmax, ρmax; step=Float64(grid_spacing))
            for y in range(-ρmax, ρmax; step=Float64(grid_spacing))
                x^2 + y^2 ≥ ρmax^2 && continue
                insidehull(SVector(x, y), hull) || continue
                pt = SVector{3,Float64}(x, y, z)
                clash = false
                for i in eachindex(atom_coords)
                    if sum(abs2, pt - atom_coords[i]) < min_dist_sq[i]
                        clash = true
                        break
                    end
                end
                clash && continue
                push!(gaussians, IsotropicGaussian(pt, Float64(σ), 1.0))
            end
        end
    end
    return IsotropicGMM(gaussians)
end
detect_pocket_cavity(seq::Chain, tmidxs; kwargs...) =
    detect_pocket_cavity(collectresidues(seq), tmidxs; kwargs...)

_complement_feature(f::Symbol) =
    f === :PosIonizable ? :NegIonizable :
    f === :NegIonizable ? :PosIonizable :
    f === :Donor        ? :Acceptor     :
    f === :Acceptor     ? :Donor        : f

"""
    mgmm = pocket_pharmacophore(seq, tmidxs; eclidxs=nothing, kwargs...)

Build a pharmacophoric `IsotropicMultiGMM` for the binding pocket of `seq`.

Inward-facing residues are identified via [`inward_tm_residues`](@ref) for the
transmembrane helices specified by `tmidxs`, and optionally via
[`inward_ecl_residues`](@ref) for the extracellular loops specified by `eclidxs`.
The pharmacophore is then built from those residues via
[`features_from_structure`](@ref); `kwargs` are forwarded to that function.

The returned `IsotropicMultiGMM` can be visualized directly with `multigmmdisplay`
from GaussianMixtureAlignment (requires a Makie backend), and the complementary
pharmacophore for ligand docking can be obtained via
[`complementary_pharmacophore`](@ref).

See also: [`complementary_pharmacophore`](@ref), [`inward_tm_residues`](@ref),
[`inward_ecl_residues`](@ref), [`features_from_structure`](@ref).
"""
function pocket_pharmacophore(seq::AbstractVector{<:AbstractResidue}, tmidxs; eclidxs=nothing, kwargs...)
    if eclidxs === nothing
        length(tmidxs) == 7 || error("Expected 7 TM helices in tmidxs, got $(length(tmidxs))")
    end
    tminward = inward_tm_residues(seq, tmidxs)
    idxs = reduce(vcat, [tm[inward] for (tm, inward) in zip(tmidxs, tminward)])
    if eclidxs === nothing
        eclidxs = default_eclidxs(tmidxs)
    end
    if !isempty(eclidxs)
        eclinward = inward_ecl_residues(seq, eclidxs)
        append!(idxs, reduce(vcat, [ecl[inward] for (ecl, inward) in zip(eclidxs, eclinward)]))
    end
    return features_from_structure(seq, idxs; kwargs...)
end
pocket_pharmacophore(seq::Chain, tmidxs; kwargs...) =
    pocket_pharmacophore(collectresidues(seq), tmidxs; kwargs...)

"""
    cmgmm = complementary_pharmacophore(mgmm)
    cmgmm = complementary_pharmacophore(mgmm, cavity; ampthreshold=0.0)

Return the pharmacophoric complement of `mgmm`, representing the features an ideal
ligand should present to bind the pocket described by `mgmm`.

Features are swapped according to the rules of molecular recognition:
- `:PosIonizable` (Arg/Lys on the receptor) ↔ `:NegIonizable` (anionic group on ligand)
- `:NegIonizable` (Asp/Glu on the receptor) ↔ `:PosIonizable` (cationic group on ligand)
- `:Donor` ↔ `:Acceptor`
- `:Steric`, `:Hydrophobe`, `:Aromatic` are unchanged

When called with a single argument, each complementary feature is placed at the same
position as the corresponding receptor atom. This is a useful approximation but the
resulting positions lie inside the receptor's van der Waals surface, so a ligand placed
there would clash with the receptor.

When `cavity` (an `IsotropicGMM` from [`detect_pocket_cavity`](@ref)) is supplied,
each complementary feature is instead placed at the nearest void-space probe position,
projecting it out of the receptor into accessible space. This form is preferred for
docking with `rocs_align` or `gogma_align` (from GaussianMixtureAlignment).

See also: [`pocket_pharmacophore`](@ref), [`detect_pocket_cavity`](@ref).

# Examples

```
julia> receptor_gmm = pocket_pharmacophore(chain, tmidxs);

julia> cavity = detect_pocket_cavity(chain, tmidxs);

julia> ligand_gmm = complementary_pharmacophore(receptor_gmm, cavity; ampthreshold=0.1);
"""
function complementary_pharmacophore(mgmm::IsotropicMultiGMM)
    result = IsotropicMultiGMM(Dict{keytype(mgmm), valtype(mgmm)}())
    for (k, gmm) in mgmm
        dest = get!(valtype(result), result, _complement_feature(k))
        for g in gmm
            push!(dest, g)
        end
    end
    return result
end

function complementary_pharmacophore(mgmm::IsotropicMultiGMM, cavity::IsotropicGMM; ampthreshold::Real=0.0)
    cavity_positions = [g.μ for g in cavity]
    result = IsotropicMultiGMM(Dict{keytype(mgmm), valtype(mgmm)}())
    for (k, gmm) in mgmm
        dest = get!(valtype(result), result, _complement_feature(k))
        for g in gmm
            nearest = argmin(sum(abs2, g.μ - c) for c in cavity_positions)
            μp = cavity_positions[nearest]
            amp = exp(-(sum(abs2, g.μ - μp) / (4 * g.σ^2)))
            push!(dest, IsotropicGaussian(cavity_positions[nearest], g.σ, amp))
        end
        # Consolidate by μ and σ before applying the amplitude threshold
        if !isempty(dest)
            g = first(dest)
            samegs = Dict{Tuple{typeof(g.μ), typeof(g.σ)}, typeof(g)}()
            for g in dest
                gagg = get(samegs, (g.μ, g.σ), nothing)
                if gagg === nothing
                    samegs[(g.μ, g.σ)] = g
                else
                    samegs[(g.μ, g.σ)] = IsotropicGaussian(g.μ, g.σ, gagg.ϕ + g.ϕ)
                end
            end
            gs = collect(Iterators.filter(g -> g.ϕ > ampthreshold, values(samegs)))
            result.gmms[k] = IsotropicGMM(gs)
        end
    end
    return result
end
