const InteractionList{K,T<:Number} = AbstractVector{Pair{Tuple{K,K},T}}
getfeatures(p::Pair{Tuple{K,K},T}) where {K,T} = p.first
getcoef(p::Pair{Tuple{K,K},T}) where {K,T} = p.second

"""
    forces = forcecomponents(seq, interactions::AbstractVector, residueindexes=eachindex(seq); kwargs...)

Calculate the forces between residues in `seq` based on the features of each
residue and the given `interactions`, which must be a list of 2-tuples
`(:field1, :field2)` or pairs `(:field1, :field2) => coef`, where `:field1` and
`:field2` are the names of the features that interact and `coef` is the
coefficient of the force (defaults to 1). The optional keyword arguments
`kwargs` are as described in [`features_from_structure`](@ref).

The feature-names can be the ones used in [`features_from_structure`](@ref).
Each pair of interactions should be listed, but only in one order (symmetry is
automatically enforced). Optionally, you can also "bundle" features together:
`:Ionic` is a bundle of `:PosIonizable` and `:NegIonizable`, where like charges
repel and opposite charges attract with the same magnitude of force. Using
`(:Ionic,:Ionic)` instead of listing all three interactions
(`:PosIonizable,:PosIonizable`, `:NegIonizable,:NegIonizable`, and
`:PosIonizable,:NegIonizable`) separately ensures that any tuning of ionic
forces will satisfy the symmetries of the real world.

Upon return, there is one force-matrix for each residue in `seq` listed in
`residueindexes`. Each force-matrix is a `3×n` matrix where each row corresponds
to a force component `(x, y, z)` and the `k`th column corresponds to
`interactions[k]`.

# Examples

```julia
julia> seq = getchain("1GZM.pdb")
interactions = [(:Steric, :Steric) => 1,           # repulsive
                (:Hydrophobe, :Hydrophobe) => -1,  # attractive
                (:Donor, :Acceptor) => -1          # attractive
                ]
julia> forces = forcecomponents(seq, interactions)
```
The output is a list of 3×4 matrices, one for each residue in `seq`.

See also: [`optimize_weights`](@ref), [`forcedict`](@ref).
"""
function forcecomponents(seq::AbstractVector{<:AbstractResidue}, interactions::InteractionList, residueindexes=eachindex(seq); kwargs...)
    # Create a MGMM for each residue
    mgmms = [IsotropicMultiGMM(Dict{Symbol,IsotropicGMM{3,Float64}}()) for _ in seq]
    for (i, seqi) in enumerate(seq)
        add_features_from_residue!(mgmms[i], seqi; kwargs...)
    end
    forces = forcecomponents(mgmms, interactions, residueindexes)
    # Project out the component of the force that would stretch or compress the covalent bonds between residues
    for (ii, i) in enumerate(residueindexes)
        force = forces[ii]
        if i > firstindex(seq)
            # Get the alpha-carbon to nitrogen vector
            v = coords(alpha_nitrogen(seq[i-1])) - coords(alpha_carbon(seq[i]))
            v /= norm(v)
            # Project out the component of the force that would stretch or compress the bond
            force -= v * (v' * force)
        end
        if i < lastindex(seq)
            # Get the nitrogen to alpha-carbon vector
            v = coords(alpha_carbon(seq[i+1])) - coords(alpha_nitrogen(seq[i]))
            v /= norm(v)
            # Project out the component of the force that would stretch or compress the bond
            force -= v * (v' * force)
        end
        forces[ii] = force
    end
    return forces
end
forcecomponents(seq, interactions::InteractionList, residueindexes=eachindex(seq); kwargs...) =
    forcecomponents(collectresidues(seq), interactions, residueindexes; kwargs...)

function forcecomponents(mgmms::AbstractVector{<:IsotropicMultiGMM{N}}, interactions::InteractionList, residueindexes=eachindex(mgmms)) where N
    interactionsdict = Dict(interactions)
    GaussianMixtureAlignment.validate_interactions(interactionsdict) || throw(ArgumentError("interactions must appear in only one order, got $interactions"))
    # We don't support (:Ionic, x) where x != :Ionic
    for interaction in interactions
        key1, key2 = getfeatures(interaction)
        if (key1 === :Ionic && key2 !== :Ionic) || (key1 !== :Ionic && key2 === :Ionic)
            throw(ArgumentError("(:Ionic, x) is not supported unless x=:Ionic, got $interaction"))
        end
    end
    cionic = get(interactionsdict, (:Ionic, :Ionic), nothing)
    if cionic !== nothing
        # If (:Ionic, :Ionic) is present, the code below will handle the signs. But we need to ensure consistency with other interactions
        for (key, sprod) in ((:Steric, :Steric) => 1, (:Hydrophobe, :Hydrophobe) => -1, (:Donor, :Acceptor) => -1)
            c = get(interactionsdict, key, nothing)
            if c !== nothing
                sign(cionic) * sign(c) != sprod && throw(ArgumentError("(:Ionic, :Ionic) is not sign-compatible with $key"))
            end
        end
    end
    forces = [zeros(N, eachindex(interactions)) for _ in residueindexes]
    for (k, interaction) in pairs(interactions)
        key1, key2 = getfeatures(interaction)
        coef = getcoef(interaction)
        # Handle feature "bundles" (e.g., PosIonizable and NegIonizable are both :Ionic)
        pairlist = key1 === key2 === :Ionic ? [(:PosIonizable, :PosIonizable) => coef, (:NegIonizable, :NegIonizable) => coef, (:PosIonizable, :NegIonizable) => -coef, (:NegIonizable, :PosIonizable) => -coef] :
                   key1 === key2 ? [(key1, key2) => coef] :
                   [(key1, key2) => coef, (key2, key1) => coef]
        for (key, c) in pairlist
            for (ii, i) in enumerate(residueindexes)
                mgmmi = mgmms[i]
                gmmx = get(mgmmi, key[1], nothing)
                gmmx === nothing && continue
                f = @view(forces[ii][:, k])
                for j in eachindex(mgmms)
                    i == j && continue
                    mgmmj = mgmms[j]
                    gmmy = get(mgmmj, key[2], nothing)
                    gmmy === nothing && continue
                    for x in gmmx, y in gmmy
                        force!(f, x, y; coef=c)
                    end
                end
            end
        end
    end
    return forces
end

"""
    w = optimize_weights(forces)

Tune the weights `w` to approximately satisfy force balance, i.e., solve

    min w.r.t. w  sum(sum(abs2, f*w) for f in forces)
    subject to    sum(w) == 1, w .>= 0

This is based on the notion that the protein structure is presumably at an
energy minimum.

This seems to work best for ubiquitous interactions, like `:Steric`, `:Hydrophobe`,
and hydrogen-bonding. Rarer interactions (`:Ionic`, `:Aromatic`) may need to be
tuned via different principles.

!!! note

    This function requires that you manually load JuMP and HiGHS, e.g., `using JuMP, HiGHS`.
"""
optimize_weights(forces::AbstractVector{<:AbstractMatrix}) = optimize_weights(sum(f' * f for f in forces))

"""
    interactiondict = forcedict(interactions::AbstractVector, w = ones(length(interactions)))

Create a dictionary of interactions, where `interactions[i]` should be weighted by `w[i]`.
"""
function forcedict(interactions::InteractionList{K}, w::AbstractVector=fill(true, eachindex(interactions))) where K
    eachindex(w) == eachindex(interactions) || throw(ArgumentError("w and interactions must have the same length"))
    # Assume `w` came from tuning the forces and is matched to `interactions`, so it shouldn't need validating
    d = Dict{Tuple{K,K},Float64}()
    for (i, interaction) in pairs(interactions)
        key1, key2 = getfeatures(interaction)
        c = getcoef(interaction)
        wi = w[i]
        if key1 === key2 === :Ionic
            d[(:PosIonizable, :PosIonizable)] =  c * wi
            d[(:NegIonizable, :NegIonizable)] =  c * wi
            d[(:PosIonizable, :NegIonizable)] = -c * wi
        else
            d[key1, key2] = c * wi
        end
    end
    return d
end
