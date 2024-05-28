const σdflt = (; hbond=1, ionic=5, aromatic=3)  # rough guesses for the interaction distance, in Å (separation is 2σ)

function default_σfun(a, r, f::Symbol; σhbond=σdflt.hbond, σionic=σdflt.ionic, σaromatic=σdflt.aromatic)
    f === :Steric && return vanderwaalsradius[(r.id.name, a.atom)]             # repulsive
    f === :Hydrophobe && return 2 * vanderwaalsradius[(r.id.name, a.atom)]     # attractive
    f ∈ (:Donor, :Acceptor) && return sqrt(σhbond^2 + vanderwaalsradius[(r.id.name, a.atom)]^2)
    f ∈ (:PosIonizable, :NegIonizable) && return σionic
    f === :Aromatic && return σaromatic # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530875/
    throw(ArgumentError("Unknown feature type: $f"))
end

equalvolumeradius(radii) = sum(radii.^3)^(1/3)
equalvolumeradius(atoms::AbstractVector{PDBAtom}, r::PDBResidue) = equalvolumeradius([vanderwaalsradius[(r.id.name, a.atom)] for a in atoms])
equalvolumeradius(atoms::AbstractVector{PDBAtom}, r::PDBResidue, ::Symbol) = equalvolumeradius(atoms, r)

function combinedfeaturesize(atoms, r, f::Symbol; σhbond=σdflt.hbond, σionic=σdflt.ionic, σaromatic=σdflt.aromatic)
    μ = mean(a -> a.coordinates, atoms)
    length(atoms) == 1 && return equalvolumeradius(atoms, r)
    σ = std([norm(a.coordinates .- μ) for a in atoms]) + equalvolumeradius(atoms, r)
    f === :Hydrophobe && return σ + 3//2
    f === :Steric && return σ
    f ∈ (:Donor, :Acceptor) && return sqrt(σhbond^2 + σ^2)
    f ∈ (:PosIonizable, :NegIonizable) && return sqrt(σ^2 + σionic^2)
    f === :Aromatic && return sqrt(σ^2 + σaromatic^2)
    throw(ArgumentError("Unknown feature type: $f"))
end

function default_ϕfun(a, r, f::Symbol)
    name = r.id.name
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

function combinedamplitude(atoms, r, f::Symbol)
    name = r.id.name
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
    if !haskey(mgmm.gmms, key)
        push!(mgmm.gmms, key => IsotropicGMM([IsotropicGaussian(vec(μ), σ, ϕ)]))
    else
        push!(mgmm.gmms[key].gaussians, IsotropicGaussian(vec(μ), σ, ϕ))
    end
end

function add_features_from_atom!(mgmm::IsotropicMultiGMM, a::PDBAtom, r::PDBResidue, σfun, ϕfun)
    ϕ = ϕfun(a, r, :Steric)
    !iszero(ϕ) && add_feature!(mgmm, :Steric, a.coordinates, σfun(a, r, :Steric), ϕ)
    if ishydrophobic(a, r.id.name)
        ϕ = ϕfun(a, r, :Hydrophobe)
        !iszero(ϕ) && add_feature!(mgmm, :Hydrophobe, a.coordinates, σfun(a, r, :Hydrophobe), ϕ)
    end
    if isaromatic(a, r.id.name)
        ϕ = ϕfun(a, r, :Aromatic)
        !iszero(ϕ) && add_feature!(mgmm, :Aromatic, a.coordinates, σfun(a, r, :Aromatic), ϕ)
    end
    if iscationic(a, r.id.name)
        ϕ = ϕfun(a, r, :PosIonizable)
        !iszero(ϕ) && add_feature!(mgmm, :PosIonizable, a.coordinates, σfun(a, r, :PosIonizable), ϕ)
    end
    if isanionic(a, r.id.name)
        ϕ = ϕfun(a, r, :NegIonizable)
        !iszero(ϕ) && add_feature!(mgmm, :NegIonizable, a.coordinates, σfun(a, r, :NegIonizable), ϕ)
    end
    if ishbonddonor(a, r.id.name)
        ϕ = ϕfun(a, r, :Donor)
        !iszero(ϕ) && add_feature!(mgmm, :Donor, a.coordinates, σfun(a, r, :Donor), ϕ)
    end
    if ishbondacceptor(a, r.id.name)
        ϕ = ϕfun(a, r, :Acceptor)
        !iszero(ϕ) && add_feature!(mgmm, :Acceptor, a.coordinates, σfun(a, r, :Acceptor), ϕ)
    end
end

function add_features_from_atoms!(mgmm::IsotropicMultiGMM, key, atoms, r, σfun, ϕfun)
    isempty(atoms) && return
    μ = mean(a -> a.coordinates, atoms)
    σ = σfun(atoms, r, key)
    ϕ = ϕfun(atoms, r, key)
    add_feature!(mgmm, key, μ, σ, ϕ)
end

function add_features_from_residue!(
    mgmm::IsotropicMultiGMM,
    r::PDBResidue;
    combined = false,
    σfun = combined ? combinedfeaturesize : default_σfun,
    ϕfun = combined ? combinedamplitude : default_ϕfun,
    weights = Dict(:Steric => 1, :Hydrophobe => 1, :Ionic => 1, :Hbond => 1, :Aromatic => 1)
)
    lookup = Dict(:PosIonizable => :Ionic, :NegIonizable => :Ionic, :Donor => :Hbond, :Acceptor => :Hbond)
    weightedϕfun(a, r, f) = ϕfun(a, r, f) * sqrt(weights[get(lookup, f, f)])   # sqrt because ϕ is squared in the overlap calculation
    if !combined
        for a in r.atoms
            add_features_from_atom!(mgmm, a, r, σfun, weightedϕfun)
        end
    else
        !iszero(weights[:Steric])     && add_features_from_atoms!(mgmm, :Steric,       r.atoms, r, σfun, weightedϕfun)
        !iszero(weights[:Hydrophobe]) && add_features_from_atoms!(mgmm, :Hydrophobe,   filter(x -> ishydrophobic(x, r.id.name),   r.atoms), r, σfun, weightedϕfun)
        !iszero(weights[:Aromatic])   && add_features_from_atoms!(mgmm, :Aromatic,     filter(x -> isaromatic(x, r.id.name),      r.atoms), r, σfun, weightedϕfun)
        !iszero(weights[:Ionic])      && add_features_from_atoms!(mgmm, :PosIonizable, filter(x -> iscationic(x, r.id.name),      r.atoms), r, σfun, weightedϕfun)
        !iszero(weights[:Ionic])      && add_features_from_atoms!(mgmm, :NegIonizable, filter(x -> isanionic(x, r.id.name),       r.atoms), r, σfun, weightedϕfun)
        !iszero(weights[:Hbond])      && add_features_from_atoms!(mgmm, :Donor,        filter(x -> ishbonddonor(x, r.id.name),    r.atoms), r, σfun, weightedϕfun)
        !iszero(weights[:Hbond])      && add_features_from_atoms!(mgmm, :Acceptor,     filter(x -> ishbondacceptor(x, r.id.name), r.atoms), r, σfun, weightedϕfun)
    end
end

"""
    features_from_structure(seq, idxs=1:length(seq); combined=false)

Construct an `IsotropicMultiGMM` from a `PDBResidue` sequence `seq` by adding features for each residue in `idxs`.

The `combined` keyword argument allows for the combination of atoms within a residue into a single feature, reducing the total number of features in the resulting model.

The `σfun` and `ϕfun` keyword arguments are functions that determine the standard deviation and amplitude of each gaussian feature, respectively.
Depending on the value of `combined`, these functions take either a single `PDBAtom` or a vector of `PDBAtom` as the first argument, and a `PDBResidue` as the second argument.
"""
function features_from_structure(seq, idxs=1:length(seq); kwargs...)
    mgmm = IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}())
    for i in idxs
        add_features_from_residue!(mgmm, seq[i]; kwargs...)
    end
    return mgmm
end
