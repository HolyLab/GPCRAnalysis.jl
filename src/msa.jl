"""
    filter_species!(msa, speciesname::AbstractString)

Remove all sequences from `msa` except those with [`species(sequencename)`](@ref) equal to `speciesname`.
"""
function filter_species!(msa, speciesname::AbstractString)
    mask = map(x -> species(x) == speciesname, sequencenames(msa))
    filtersequences!(msa, mask)
end

"""
    filter_long!(msa, minres::Real)

Remove all sequences from `msa` with fewer than `minres` matching residues.
"""
function filter_long!(msa, minres::Real)
    # Get rid of short sequences
    nresidues = map(eachrow(msa)) do v
        sum(!=(Residue('-')), v)
    end
    mask = nresidues .> minres
    filtersequences!(msa, mask)
end
