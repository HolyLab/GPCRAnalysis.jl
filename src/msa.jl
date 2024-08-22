Base.getindex(msa::AbstractMultipleSequenceAlignment, seqname::MSACode) = getsequence(msa, seqname.name)
Base.getindex(msa::AbstractMultipleSequenceAlignment, seqname::AccessionCode) = getsequence(msa, MSACode(msa, seqname).name)

"""
    tour = sortperm_msa(msa::AbstractMultipleSequenceAlignment)

Order the sequences in `msa` to minimize the "tour length" visiting each sequence once.
The length between sequences is defined at `100 - percentsimilarity(seq1, seq2)`.

This can be useful for graphical or alignment display by grouping obviously-similar
sequences near one another.
"""
function sortperm_msa(msa::AbstractMultipleSequenceAlignment)
    sim = percentsimilarity(msa)
    D = 100 .- Matrix(sim)
    tour, _ = solve_tsp(D; quality_factor=80)
    pop!(tour)  # returns a closed path, but that duplicates the first/last entry
    return tour
end

"""
    filter_species!(msa::AbstractMultipleSequenceAlignment, speciesname::AbstractString)

Remove all sequences from `msa` except those with [`species(sequencename)`](@ref) equal to `speciesname`.
"""
function filter_species!(msa::AbstractMultipleSequenceAlignment, speciesname::AbstractString)
    mask = map(x -> species(x) == speciesname, sequencenames(msa))
    filtersequences!(msa, mask)
end

"""
    filter_long!(msa::AbstractMultipleSequenceAlignment, minres::Real)

Remove all sequences from `msa` with fewer than `minres` matching residues.
"""
function filter_long!(msa::AbstractMultipleSequenceAlignment, minres::Real)
    # Get rid of short sequences
    nresidues = map(eachrow(msa)) do v
        sum(!=(MSA.Residue('-')), v)
    end
    mask = nresidues .> minres
    filtersequences!(msa, mask)
end

"""
    ac = AccessionCode(msa, seqname)

Return the Uniprot accession code associated with `seqname`.
"""
function AccessionCode(msa::AnnotatedMultipleSequenceAlignment, seqname::AbstractString)
    AccessionCode(uniprotX(getannotsequence(msa, seqname, "AC", seqname)))
end
AccessionCode(msa::AnnotatedMultipleSequenceAlignment, seqname::MSACode) = AccessionCode(msa, seqname.name)
AccessionCode(::AnnotatedMultipleSequenceAlignment, seqname::AccessionCode) = seqname

function MSACode(msa::AnnotatedMultipleSequenceAlignment, accession::AbstractString)
    seqnames = sequencenames(msa)
    return MSACode(seqnames[findfirst(x -> AccessionCode(msa, x).name == accession, seqnames)])
end
MSACode(msa::AnnotatedMultipleSequenceAlignment, accession::AccessionCode) = MSACode(msa, accession.name)
MSACode(::AnnotatedMultipleSequenceAlignment, accession::MSACode) = accession

# Move this to MIToS?
if !hasmethod(getsequencemapping, Tuple{AnnotatedAlignedSequence})
    function MIToS.MSA.getsequencemapping(seq::AnnotatedAlignedSequence)
        getsequencemapping(seq, sequencenames(seq)[1])
    end
    function MIToS.MSA.getsequencemapping(msa::Union{AnnotatedAlignedSequence,AnnotatedMultipleSequenceAlignment}, seq_id::String)
        MIToS.MSA._str2int_mapping(getannotsequence(msa, seq_id, "SeqMap"))
    end
    function MIToS.MSA.getsequencemapping(msa::AnnotatedMultipleSequenceAlignment, seqid::Regex)
        id = findfirst(str -> occursin(seqid, str), sequencenames(msa))
        getsequencemapping(msa, id)
    end
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
SequenceMapping(seq::AnnotatedAlignedSequence) = SequenceMapping(getsequencemapping(seq))

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
