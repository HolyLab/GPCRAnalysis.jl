const σdflt = (; hbond=1, ionic=5, aromatic=3)  # rough guesses for the interaction distance, in Å (separation is 2σ)

function stdname(aname)
    if aname == "OXT"
        return "O"
    end
    return aname
end

equalvolumeradius(radii) = sum(radii.^3)^(1/3)
equalvolumeradius(atoms::ResidueLike, r::AbstractResidue) = equalvolumeradius([vanderwaalsradius[(resname(r), stdname(atomname(a)))] for a in atoms])
equalvolumeradius(atoms::ResidueLike, r::AbstractResidue, ::Symbol) = equalvolumeradius(atoms, r)

function atomic_σfun(a, r, f::Symbol; σhbond=σdflt.hbond, σionic=σdflt.ionic, σaromatic=σdflt.aromatic)
    aname, rname = stdname(atomname(a)), resname(r)
    f === :Steric && return vanderwaalsradius[(rname, aname)]             # repulsive
    f === :Hydrophobe && return 2 * vanderwaalsradius[(rname, aname)]     # attractive
    f ∈ (:Donor, :Acceptor) && return sqrt(σhbond^2 + vanderwaalsradius[(rname, aname)]^2)
    f ∈ (:PosIonizable, :NegIonizable) && return σionic
    f === :Aromatic && return σaromatic # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530875/
    throw(ArgumentError("Unknown feature type: $f"))
end

function combined_σfun(atoms, r, f::Symbol; σhbond=σdflt.hbond, σionic=σdflt.ionic, σaromatic=σdflt.aromatic)
    μ = mean(a -> coords(a), atoms)
    length(atoms) == 1 && return equalvolumeradius(atoms, r)
    σ = std([norm(coords(a) .- μ) for a in atoms]) + equalvolumeradius(atoms, r)
    f === :Hydrophobe && return σ + 3//2
    f === :Steric && return σ
    f ∈ (:Donor, :Acceptor) && return sqrt(σhbond^2 + σ^2)
    f ∈ (:PosIonizable, :NegIonizable) && return sqrt(σ^2 + σionic^2)
    f === :Aromatic && return sqrt(σ^2 + σaromatic^2)
    throw(ArgumentError("Unknown feature type: $f"))
end

function atomic_ϕfun(a, r, f::Symbol)
    name = resname(r)
    if f === :PosIonizable
        return aggcharges[name][atomname(a)]
    end
    if f === :NegIonizable
        return abs(aggcharges[name][atomname(a)])
    end
    return 1.0
end

function combined_ϕfun(atoms, r, f::Symbol)
    name = resname(r)
    if f === :PosIonizable
        name == "HIS" && return 0.1
        name == "ARG" && return 1.0
        name == "LYS" && return 1.0
        return 0.0
    end
    if f === :NegIonizable
        name ∈ ("ASP", "GLU") && return 1.0
        return 0.0
    end
    f === :Aromatic && return 1.0
    return Float64(length(atoms))
end

function add_feature!(mgmm::IsotropicMultiGMM, key, μ, σ, ϕ)
    gmm = get!(valtype(mgmm), mgmm, key)
    push!(gmm, eltype(gmm)(vec(μ), σ, ϕ))
end

function add_features_from_atom!(mgmm::IsotropicMultiGMM, a::AbstractAtom, r::AbstractResidue, σfun, ϕfun)
    ϕ = ϕfun(a, r, :Steric)
    rname, aname, acoords = resname(r), atomname(a), coords(a)
    first(aname) == 'H' && return mgmm # skip hydrogens; their properties are included in the heavy atom they are bonded to
    !iszero(ϕ) && add_feature!(mgmm, :Steric, acoords, σfun(a, r, :Steric), ϕ)
    z = if rname == "HIS"
        return aname == "NE2" ? 0.1 : 0.0
    else
        get(aggcharges[rname], aname, 0.0f0)
    end
    if abs(z) < 0.3
        ϕ = ϕfun(a, r, :Hydrophobe)
        !iszero(ϕ) && add_feature!(mgmm, :Hydrophobe, acoords, σfun(a, r, :Hydrophobe), ϕ)
    end
    if isaromatic(aname, rname)
        ϕ = ϕfun(a, r, :Aromatic)
        !iszero(ϕ) && add_feature!(mgmm, :Aromatic, acoords, σfun(a, r, :Aromatic), ϕ)
    end
    if z >= 0.3
        ϕ = z
        !iszero(ϕ) && add_feature!(mgmm, :PosIonizable, acoords, σfun(a, r, :PosIonizable), ϕ)
    end
    if z <= -0.3
        ϕ = z
        !iszero(ϕ) && add_feature!(mgmm, :NegIonizable, acoords, σfun(a, r, :NegIonizable), ϕ)
    end
    if ishbonddonor(aname, rname)
        ϕ = ϕfun(a, r, :Donor)
        !iszero(ϕ) && add_feature!(mgmm, :Donor, acoords, σfun(a, r, :Donor), ϕ)
    end
    if ishbondacceptor(aname, rname)
        ϕ = ϕfun(a, r, :Acceptor)
        !iszero(ϕ) && add_feature!(mgmm, :Acceptor, acoords, σfun(a, r, :Acceptor), ϕ)
    end
    return mgmm
end

function add_features_from_atoms!(mgmm::IsotropicMultiGMM, key, atoms, r, σfun, ϕfun)
    isempty(atoms) && return
    μ = mean(coords, atoms)
    σ = σfun(atoms, r, key)
    ϕ = ϕfun(atoms, r, key)
    add_feature!(mgmm, key, μ, σ, ϕ)
end

function add_features_from_residue!(
    mgmm::IsotropicMultiGMM,
    r::AbstractResidue;
    combined = false,
    σfun = combined ? combined_σfun : atomic_σfun,
    ϕfun = combined ? combined_ϕfun : atomic_ϕfun
)
    if !combined
        for a in r
            add_features_from_atom!(mgmm, a, r, σfun, ϕfun)
        end
    else
        rname = resname(r)
        add_features_from_atoms!(mgmm, :Steric,       r, r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Hydrophobe,   collectatoms(r, x -> ishydrophobic(atomname(x), rname))  , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Aromatic,     collectatoms(r, x -> isaromatic(atomname(x), rname))     , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :PosIonizable, collectatoms(r, x -> isstdcationic(atomname(x), rname))     , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :NegIonizable, collectatoms(r, x -> isstdanionic(atomname(x), rname))      , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Donor,        collectatoms(r, x -> ishbonddonor(atomname(x), rname))   , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Acceptor,     collectatoms(r, x -> ishbondacceptor(atomname(x), rname)), r, σfun, ϕfun)
    end
end

"""
    mgmm = features_from_structure(seq::ChainLike, idxs=1:length(seq); combined=false)

Construct an `IsotropicMultiGMM` from `seq` by adding
features for each residue in `idxs`.

`combined=true` causes all atoms sharing the same feature to be combined,
reducing the total number of features in the resulting model. The default
creates features for each atom separately.

The `σfun` and `ϕfun` keyword arguments are functions that determine the
standard deviation and amplitude of each gaussian feature, respectively,
and take arguments `(atom, residue, feature)`.

The output `mgmm` may include the following features:
- `:Steric` (the "hard center", a proxy for the repulsive core of Lennard-Jones potentials)
- `:Hydrophobe` (van der Waals interactions, a proxy for the attractive part of Lennard-Jones potentials)
- `:Aromatic` (the aromatic ring of phenylalanine, tyrosine, and tryptophan)
- `:PosIonizable` (the positively charged nitrogen of histidine, arginine, and lysine; histidine gets a fractional charge of +0.1)
- `:NegIonizable` (the negatively charged oxygen of aspartate and glutamate; each oxygen gets a fractional charge of -0.5)
- `:Donor` (the hydrogen of a hydrogen bond donor)
- `:Acceptor` (the oxygen of a hydrogen bond acceptor)
"""
function features_from_structure(seq, idxs=1:length(seq); kwargs...)
    mgmm = IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}())
    for i in idxs
        add_features_from_residue!(mgmm, seq[i]; kwargs...)
    end
    return mgmm
end

"""
    mgmm = features_from_structure(seq::ChainLike, ρmax::Real, zi::AbstractInterval)

Construct an `IsotropicMultiGMM` from `seq` including all atoms that lie within
the cylinder

    x^2 + y^2 <= ρmax^2
    z ∈ zi

`zi` is an `AbstractInterval`, e.g. `0..30` or `-15..15` (see IntervalSets.jl).

This implicitly assumes that you've aligned `seq` to the membrane, or aligned
`seq` to a homolog that is membrane-aligned. See [`align_to_membrane`](@ref),
[`align`](@ref).
"""
function features_from_structure(
        seq, ρmax::Real, zi::AbstractInterval;
        σfun = atomic_σfun,
        ϕfun = atomic_ϕfun,
    )
    mgmm = IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}())
    for r in seq
        for a in r
            c = coords(a)
            ρ = hypot(c[1], c[2])
            if ρ < ρmax && c[3] ∈ zi
                add_features_from_atom!(mgmm, a, r, σfun, ϕfun)
            end
        end
    end
    return mgmm
end
