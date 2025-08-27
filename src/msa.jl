## Low level MSA API
## Different MSA extensions (MIToS or BioStockholm) should implement these

# It's required that `Char(res)` converts a residue to a `Char`

"""
    idxs = sequenceindexes(msaseq)
    idxs = sequenceindexes(msa, i)

Return the corresponding index within the full sequence for each position in `msaseq`.
`0` indicates a gap or unknown residue.

The two-argument form retrieves the sequenceindexes for the `i`th sequence in `msa`.
"""
function sequenceindexes end
sequenceindexes(msa::AbstractVector{FASTX.FASTA.Record}, i::Int) = sequenceindexes(msa[i])
function sequenceindexes(seq::FASTX.FASTA.Record)
    n = 0
    return [isgap(c) ? 0 : (n += 1) for c in sequence(seq)]
end

"""
    idxs = columnindexes(msa)

Return the indices of the conserved columns of the MSA.
"""
function columnindexes end

"""
    isgap(res)

Return `true` if the residue `res` is a gap.
"""
function isgap end
isgap(c::Char) = c == '-'
isgap(r::Integer) = iszero(r)   # when mapped to integers, gaps are encoded as 0

"""
    isunknown(res)

Return `true` if the residue `res` is unknown.
"""
function isunknown end
isunknown(c::Char) = c == 'X'

"""
    keys = sequencekeys(msa)

Return the keys (sequence names) of the MSA.
"""
function sequencekeys end
sequencekeys(msa::AbstractVector{FASTX.FASTA.Record}) = eachindex(msa)

"""
    seq = msasequence(msa, key)

Return the aligned sequence corresponding to `key`.
"""
function msasequence end
msasequence(msa::AbstractVector{FASTX.FASTA.Record}, key::Int) = sequence(msa[key])

"""
    R = residuematrix(msa)

Get all residues in the MSA as a matrix, one sequence per row.
"""
function residuematrix end
function residuematrix(msa::AbstractVector{FASTX.FASTA.Record})
    reduce(vcat, [permutedims(collect(sequence(rec))) for rec in msa])
end

"""
    msaview = subseqs(msa, rowindexes::AbstractVector{Int})
    msaview = subseqs(msa, rowmask::AbstractVector{Bool})

    subseqs!(msa, rowindexes::AbstractVector{Int})
    subseqs!(msa, rowmask::AbstractVector{Bool})

Construct a reduced-size `msaview`, keeping only the sequences corresponding to `rowindexes`/`rowmask`.
"""
function subseqs end
function subseqs! end

## End required API, but some can specialize other methods

"""
    pc = percent_similarity(msa)
    pc = percent_similarity(f, msa)

Compute the percent similarity between all pairs of sequences in `msa`.
`pc[i, j]` is the percent similarity between sequences `i` and `j`.

Optionally apply mapping function `f` to each residue before computing
similarity.
"""
function percent_similarity end

function percent_similarity(f, msa)
    # This mimics MIToS's implementation
    function pctsim(v1, v2)
        same = l = 0
        for (a, b) in zip(v1, v2)
            isgap(a) && isgap(b) && continue  # skip gaps
            same += a == b
            l += 1
        end
        return 100 * same / l
    end

    M = f.(residuematrix(msa))
    n = size(M, 1)
    S = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            S[i, j] = pctsim(M[i, :], M[j, :])
            S[j, i] = S[i, j]
        end
    end
    return S
end
percent_similarity(msa) = percent_similarity(reduced_alphabet, msa)

function reduced_alphabet(r::Char)
    if r == '-'
        return 0
    elseif r in ('A','I','L','M','V')
        return 1  # hydrophobic
    elseif r in ('N','Q','S','T')
        return 2  # polar
    elseif r in ('R','H','K')
        return 3  # charged
    elseif r in ('D','E')
        return 4  # charged
    elseif r in ('F','W','Y')
        return 5  # aromatic
    end
    offset = findfirst(==(r), ('C','G','P'))
    offset === nothing && throw(ArgumentError("Unknown residue '$r'"))
    return 5 + offset  # special or unknown
end

columnwise_entropy(msa) = columnwise_entropy(reduced_alphabet, msa)

# Notes on interpreting letter codes in the "GC.seq_cons" field:
# - `.` indicates a gap
# - uppercase single-letter amino acid codes indicate strong consensus (>60%)
# - lowercase single-letter codes (likely interpretations):
#   + 'a': aromatic
#   + 'h': hydrophobic
#   + rest unknown
# - '+' and '-' indicate positively- and negatively-charged residues, respectively

## MSA functions

"""
    tour = sortperm_msa(msa)

Order the sequences in `msa` to minimize the "tour length" visiting each sequence once.
The length between sequences is defined at `100 - percentsimilarity(seq1, seq2)`.

This can be useful for graphical or alignment display by grouping obviously-similar
sequences near one another.
"""
function sortperm_msa(msa)
    sim = percent_similarity(msa)
    D = 100 .- Matrix(sim)
    tour, _ = solve_tsp(D; quality_factor=80)
    pop!(tour)  # returns a closed path, but that duplicates the first/last entry
    return tour
end

"""
    filter_species!(msa, speciesname::AbstractString)

Remove all sequences from `msa` except those with [`species(sequencename)`](@ref) equal to `speciesname`.
"""
function filter_species!(msa, speciesname::AbstractString)
    mask = map(x -> species(x) == speciesname, sequencekeys(msa))
    subseqs!(msa, mask)
end

"""
    filter_long!(msa, minres::Real)

Remove all sequences from `msa` with fewer than `minres` matching residues.
"""
function filter_long!(msa, minres::Real)
    # Get rid of short sequences
    nresidues = map(sequencekeys(msa)) do key
        sum(!isgap, msasequence(msa, key))
    end
    rowmask = nresidues .> minres
    subseqs!(msa, rowmask)
end

struct SequenceMapping <: AbstractVector{Int}
    seqmap::Vector{Int}
end

"""
    sm = SequenceMapping([4, 5, 0, ...])
    sm = SequenceMapping(seq::AnnotatedAlignedSequence)

A `SequenceMapping` is a vector of indexes within a full sequence that map to a
reference. Specifically, `sm[i]` is the index of the residue in the full
sequence that maps to the `i`-th position in the reference. 0 is a placeholder
for a position in the reference that has no mapping to the full sequence.

# Example

`SequenceMapping([4, 5, 0, ...])` indicates that:
- the first position in the reference maps to the fourth residue in the full sequence,
- the second position in the reference maps to the fifth residue in the full sequence, and
- the third position in the reference lacks a corresponding residue in the full sequence.
"""
function SequenceMapping end

Base.size(sm::SequenceMapping) = size(sm.seqmap)
Base.getindex(sm::SequenceMapping, i::Int) = sm.seqmap[i]
Base.setindex!(sm::SequenceMapping, v, i::Int) = sm.seqmap[i] = v
Base.similar(::SequenceMapping, ::Type{Int}, dims::Dims{1}) = SequenceMapping(Vector{Int}(undef, dims))

# We restrict the following to `Vector` to prevent ambiguities
Base.getindex(fullseqvec::Vector{T}, sm::SequenceMapping) where T =
    Union{T,Nothing}[i == 0 ? nothing : fullseqvec[i] for i in sm.seqmap]

Base.getindex(fullseq::Chain, sm::SequenceMapping) = getindex(collectresidues(fullseq), sm)
Base.getindex(fullseq::Model, sm::SequenceMapping) = getindex(only(fullseq), sm)
Base.getindex(fullseq::MolecularStructure, sm::SequenceMapping) = getindex(only(fullseq), sm)

function Base.setindex!(fullseqvec::Vector, colvalues::AbstractVector, sm::SequenceMapping)
    eachindex(colvalues) === eachindex(sm) || throw(DimensionMismatch("`colvalues` and `sm` must have the same indices"))
    for (v, i) in zip(colvalues, sm.seqmap)
        i == 0 && continue
        fullseqvec[i] = v
    end
    return fullseqvec
end
