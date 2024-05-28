"""
    forces = forcecomponents(seq::Vector{PDBResidue}, keep::AbstractVector{Bool}; kwargs...)

Calculate the forces between residues in a `PDBResidue` sequence `seq` based on
the features of each residue. `kwargs` are passed to
`add_features_from_residue!`.

Upon return, there is one force-matrix for each residue in `seq` for which
`keep` is `true`. Each force-matrix is a 3x5 matrix where each row corresponds
to a force component (x, y, z) and each column corresponds to a feature type
(steric, hydrophobic, ionic, hbond, aromatic).
"""
function forcecomponents(seq::Vector{PDBResidue}, keep::AbstractVector{Bool}; kwargs...)
    # Create a MGMM for each residue
    mgmms = [IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}()) for _ in seq]
    for i in eachindex(seq)
        add_features_from_residue!(mgmms[i], seq[i]; kwargs...)
    end
    # steric = 1, hydrophobic = 2, ionic = 3, hbond = 4, aromatic = 5
    forces = Matrix{Float64}[]
    for i in eachindex(seq)
        keep[i] || continue
        force = zeros(3, 5)
        mgmmi = mgmms[i]
        for (colidx, keyi) in ((1, :Steric), (2, :Hydrophobe), (3, :PosIonizable), (3, :NegIonizable), (4, :Donor), (4, :Acceptor), (5, :Aromatic))
            haskey(mgmmi.gmms, keyi) || continue
            gx = mgmmi.gmms[keyi]
            for j in eachindex(seq)
                if i == j
                    continue
                end
                mgmmj = mgmms[j]
                if keyi ∈ (:Steric, :Hydrophobe, :Aromatic)
                    keyj = keyi
                    haskey(mgmmj.gmms, keyj) || continue
                    force!(@view(force[:, colidx]), gx, mgmmj.gmms[keyj]; coef = keyi === :Steric ? -1 : 1)
                elseif keyi ∈ (:PosIonizable, :NegIonizable)
                    for keyj in (:PosIonizable, :NegIonizable)
                        haskey(mgmmj.gmms, keyj) || continue
                        force!(@view(force[:, colidx]), gx, mgmmj.gmms[keyj]; coef = (keyj == keyi ? -1 : 1))
                    end
                elseif keyi === :Donor
                    force!(@view(force[:, colidx]), gx, mgmmj.gmms[:Acceptor]; coef = 1)
                elseif keyi === :Acceptor
                    force!(@view(force[:, colidx]), gx, mgmmj.gmms[:Donor]; coef = 1)
                end
            end
        end
        # Project out the component of the force that would stretch or compress the covalent bonds between residues
        if i > 1
            # Get the alpha-carbon to nitrogen vector
            v = alpha_nitrogen(seq[i-1]).coordinates - alpha_carbon(seq[i]).coordinates
            v /= norm(v)
            # Project out the component of the force that would stretch or compress the bond
            force -= v * (v' * force)
        end
        if i < length(seq)
            # Get the nitrogen to alpha-carbon vector
            v = alpha_carbon(seq[i+1]).coordinates - alpha_nitrogen(seq[i]).coordinates
            v /= norm(v)
            # Project out the component of the force that would stretch or compress the bond
            force -= v * (v' * force)
        end
        push!(forces, force)
    end
    return forces
end

function force_vector2dict(w::AbstractVector)
    length(w) == 5 && return Dict(:Steric => w[1], :Hydrophobe => w[2], :Ionic => w[3], :Hbond => w[4], :Aromatic => w[5])
    length(w) == 3 && return Dict(:Steric => w[1], :Hydrophobe => w[2], :Ionic => 0.0, :Hbond => w[3], :Aromatic => 0.0)
    error("no standard conversion from vector to dict for length $(length(w))")
end
force_dict2vector(d::Dict) = [d[:Steric], d[:Hydrophobe], d[:Ionic], d[:Hbond], d[:Aromatic]]
