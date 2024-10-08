const σdflt = (; hbond=1, ionic=5, aromatic=3)  # rough guesses for the interaction distance, in Å (separation is 2σ)

equalvolumeradius(radii) = sum(radii.^3)^(1/3)
equalvolumeradius(atoms::ResidueLike, r::AbstractResidue) = equalvolumeradius([vanderwaalsradius[(resname(r), atomname(a))] for a in atoms])
equalvolumeradius(atoms::ResidueLike, r::AbstractResidue, ::Symbol) = equalvolumeradius(atoms, r)

function atomic_σfun(a, r, f::Symbol; σhbond=σdflt.hbond, σionic=σdflt.ionic, σaromatic=σdflt.aromatic)
    aname, rname = atomname(a), resname(r)
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
        name == "HIS" && return 0.1
        name == "ARG" && return 0.5  # arginine has two ionizable atoms, but the net charge is 1
        name == "LYS" && return 1.0  # lysine has one ionizable atom
    end
    if f === :NegIonizable
        name ∈ ("ASP", "GLU") && return 0.5  # aspartate and glutamate have two
    end
    return 1.0
end

function combined_ϕfun(atoms, r, f::Symbol)
    name = resname(r)
    if f === :PosIonizable
        name == "HIS" && return 0.1
        name == "ARG" && return 1.0
        name == "LYS" && return 1.0
        throw(ArgumentError("Unknown residue type: $name"))
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
    rname, acoords = resname(r), coords(a)
    !iszero(ϕ) && add_feature!(mgmm, :Steric, acoords, σfun(a, r, :Steric), ϕ)
    if ishydrophobic(a, rname)
        ϕ = ϕfun(a, r, :Hydrophobe)
        !iszero(ϕ) && add_feature!(mgmm, :Hydrophobe, acoords, σfun(a, r, :Hydrophobe), ϕ)
    end
    if isaromatic(a, rname)
        ϕ = ϕfun(a, r, :Aromatic)
        !iszero(ϕ) && add_feature!(mgmm, :Aromatic, acoords, σfun(a, r, :Aromatic), ϕ)
    end
    if iscationic(a, rname)
        ϕ = ϕfun(a, r, :PosIonizable)
        !iszero(ϕ) && add_feature!(mgmm, :PosIonizable, acoords, σfun(a, r, :PosIonizable), ϕ)
    end
    if isanionic(a, rname)
        ϕ = ϕfun(a, r, :NegIonizable)
        !iszero(ϕ) && add_feature!(mgmm, :NegIonizable, acoords, σfun(a, r, :NegIonizable), ϕ)
    end
    if ishbonddonor(a, rname)
        ϕ = ϕfun(a, r, :Donor)
        !iszero(ϕ) && add_feature!(mgmm, :Donor, acoords, σfun(a, r, :Donor), ϕ)
    end
    if ishbondacceptor(a, rname)
        ϕ = ϕfun(a, r, :Acceptor)
        !iszero(ϕ) && add_feature!(mgmm, :Acceptor, acoords, σfun(a, r, :Acceptor), ϕ)
    end
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
        add_features_from_atoms!(mgmm, :Hydrophobe,   collectatoms(r, x -> ishydrophobic(x, rname))  , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Aromatic,     collectatoms(r, x -> isaromatic(x, rname))     , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :PosIonizable, collectatoms(r, x -> iscationic(x, rname))     , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :NegIonizable, collectatoms(r, x -> isanionic(x, rname))      , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Donor,        collectatoms(r, x -> ishbonddonor(x, rname))   , r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Acceptor,     collectatoms(r, x -> ishbondacceptor(x, rname)), r, σfun, ϕfun)
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
