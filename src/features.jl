equalvolumeradius(radii) = sum(radii.^3)^(1/3)

function rescast(f, res::PDBResidue)
    return [f(a,res.id.name) for a in res.atoms]
end

function add_feature!(mgmm::IsotropicMultiGMM, key, μ, σ, ϕ)
    if !haskey(mgmm.gmms, key)
        push!(mgmm.gmms, key => IsotropicGMM([IsotropicGaussian(μ, σ, ϕ)]))
    else
        push!(mgmm.gmms[key].gaussians, IsotropicGaussian(μ, σ, ϕ))
    end
end

function add_features_from_atom!(mgmm::IsotropicMultiGMM, a::PDBAtom, r::PDBResidue, ϕ=1.0)
    σ = vanderwaalsradius[(r.id.name, a.atom)]
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

function add_features_from_atoms!(mgmm::IsotropicMultiGMM, key, atoms, r, ϕ=1.0)
    if isempty(atoms)
        return
    end
    μ = mean(a -> a.coordinates, atoms)
    σ = equalvolumeradius([vanderwaalsradius[(r.id.name, a.atom)] for a in atoms])
    add_feature!(mgmm, key, μ, σ, ϕ)
end

function add_features_from_residue!(mgmm::IsotropicMultiGMM, r::PDBResidue, ϕ=1.0; combined = false)
    if !combined
        for a in r.atoms
            add_features_from_atom!(mgmm, a, r, ϕ)
        end
    else
        add_features_from_atoms!(mgmm, :Hydrophobe, r.atoms[rescast(ishydrophobic,r)], r, ϕ)
        add_features_from_atoms!(mgmm, :Aromatic, r.atoms[rescast(isaromatic,r)], r, ϕ)
        add_features_from_atoms!(mgmm, :PosIonizable, r.atoms[rescast(iscationic,r)], r, ϕ)
        add_features_from_atoms!(mgmm, :NegIonizable, r.atoms[rescast(isanionic,r)], r, ϕ)
        add_features_from_atoms!(mgmm, :Donor, r.atoms[rescast(ishbonddonor,r)], r, ϕ)
        add_features_from_atoms!(mgmm, :Acceptor, r.atoms[rescast(ishbondacceptor,r)], r, ϕ)
    end
end

"""
    features_from_structure(seq, idxs=1:length(seq); combined=false)

Construct an `IsotropicMultiGMM` from a `PDBResidue` sequence `seq` by adding features for each residue in `idxs`.

The `combined` keyword argument allows for the combination of atoms within a residue into a single feature, reducing the total number of features in the resulting model.
"""
function features_from_structure(seq, idxs=1:length(seq); kwargs...)
    mgmm = IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}())
    for i in idxs
        add_features_from_residue!(mgmm, seq[i]; kwargs...)
    end
    return mgmm 
end