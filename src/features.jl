equalvolumeradius(radii) = sum(radii.^3)^(1/3)
equalvolumeradius(atoms::AbstractVector{PDBAtom}, r::PDBResidue) = equalvolumeradius([vanderwaalsradius[(r.id.name, a.atom)] for a in atoms])

function combinedfeaturesize(atoms, r)
    μ = mean(a -> a.coordinates, atoms)
    length(atoms) == 1 && return equalvolumeradius(atoms, r)
    return std([norm(a.coordinates .- μ) for a in atoms]) + equalvolumeradius(atoms, r)
end

function add_feature!(mgmm::IsotropicMultiGMM, key, μ, σ, ϕ)
    if !haskey(mgmm.gmms, key)
        push!(mgmm.gmms, key => IsotropicGMM([IsotropicGaussian(vec(μ), σ, ϕ)]))
    else
        push!(mgmm.gmms[key].gaussians, IsotropicGaussian(vec(μ), σ, ϕ))
    end
end

function add_features_from_atom!(mgmm::IsotropicMultiGMM, a::PDBAtom, r::PDBResidue, σfun, ϕfun)
    σ = σfun(a, r)
    ϕ = ϕfun(a, r)
    if ishydrophobic(a, r.id.name)
        add_feature!(mgmm, :Hydrophobe, a.coordinates, σ, ϕ)
    end
    if isaromatic(a, r.id.name)
        add_feature!(mgmm, :Aromatic, a.coordinates, σ, ϕ)
    end
    if iscationic(a, r.id.name)
        add_feature!(mgmm, :PosIonizable, a.coordinates, σ, ϕ)
    end
    if isanionic(a, r.id.name)
        add_feature!(mgmm, :NegIonizable, a.coordinates, σ, ϕ)
    end
    if ishbonddonor(a, r.id.name)
        add_feature!(mgmm, :Donor, a.coordinates, σ, ϕ)
    end
    if ishbondacceptor(a, r.id.name)
        add_feature!(mgmm, :Acceptor, a.coordinates, σ, ϕ)
    end
end

function add_features_from_atoms!(mgmm::IsotropicMultiGMM, key, atoms, r, σfun, ϕfun)
    isempty(atoms) && return
    μ = mean(a -> a.coordinates, atoms)
    σ = σfun(atoms, r)
    ϕ = ϕfun(atoms, r)
    add_feature!(mgmm, key, μ, σ, ϕ)
end

function add_features_from_residue!(
    mgmm::IsotropicMultiGMM, 
    r::PDBResidue; 
    combined = false, 
    σfun = combined ? combinedfeaturesize : (a, r) -> vanderwaalsradius[(r.id.name, a.atom)],
    ϕfun = combined ? (atoms, r) -> length(atoms) : (args...) -> 1.0
)
    if !combined
        for a in r.atoms
            add_features_from_atom!(mgmm, a, r, σfun, ϕfun)
        end
    else
        add_features_from_atoms!(mgmm, :Hydrophobe,   filter(x -> ishydrophobic(x, r.id.name),   r.atoms), r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Aromatic,     filter(x -> isaromatic(x, r.id.name),      r.atoms), r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :PosIonizable, filter(x -> iscationic(x, r.id.name),      r.atoms), r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :NegIonizable, filter(x -> isanionic(x, r.id.name),       r.atoms), r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Donor,        filter(x -> ishbonddonor(x, r.id.name),    r.atoms), r, σfun, ϕfun)
        add_features_from_atoms!(mgmm, :Acceptor,     filter(x -> ishbondacceptor(x, r.id.name), r.atoms), r, σfun, ϕfun)
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

function aa_features_from_structure(seq, idxs=1:length(seq); kwargs...)
    mgmm = IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}())
    for i in idxs
        add_features_from_residue!(mgmm, seq[i]; kwargs...)
    end
    return mgmm 
end