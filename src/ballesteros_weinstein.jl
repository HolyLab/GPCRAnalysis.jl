struct BWScheme
    conserved::NTuple{7,Int}
    tmspans::NTuple{7,UnitRange{Int}}
end

"""
    BWScheme(conserved_idx, tmspans)

Specify the Ballesteros-Weinstein scheme used for a particular protein.
`conserved_idx` is a list of 7 "most conserved" residues per helix
(rhodopsin family: N1, D2, R3, W4, P5, P6, P7) and the span of each helix.

# Examples

For mouse rhodopsin (P15409),

```julia
julia> opsd_scheme = BWScheme([55, 83, 135, 161, 215, 267, 303],
                              [37:61, 74:96, 111:133, 153:173, 203:224, 253:274, 287:308]);
```
"""
function BWScheme(bwconserved::AbstractVector, tmspans::AbstractVector)
    length(bwconserved) == length(tmspans) == 7 || error("GPCRs have 7 transmembrane regions")
    return BWScheme((bwconserved...,), (tmspans...,))
end

"""
    helix, residue_position = lookupbw(idx::Integer, scheme::BWScheme)

Calculate the Ballesteros-Weinstein residue number coresponding to residue `idx`.
`tmspans` describes the transmembrane regions in the reference, and `bwconserved`
the index of the most-conserved residue.

# Examples

For mouse rhodopsin (P15409),

```julia
julia> opsd_scheme = BWScheme([55, 83, 135, 161, 215, 267, 303],
           [34:64, 73:99, 107:139, 150:173, 200:229, 246:277, 285:309]);

julia> lookupbw(160, opsd_scheme)
(4, 49)

julia> lookupbw((4, 49), opsd_scheme)
160
```

This is the residue just before the most-conserved residue of helix 4.
"""
function lookupbw(idx::Integer, scheme::BWScheme)
    i = findfirst(rng -> idx ∈ rng, scheme.tmspans)
    i === nothing && return nothing
    @assert findnext(rng -> idx ∈ rng, scheme.tmspans, i+1) === nothing
    return i, idx - scheme.conserved[i] + 50
end

function lookupbw((helix, offset)::Tuple{Integer,Integer}, scheme::BWScheme)
    icons = scheme.conserved[helix]
    return icons - 50 + offset
end
