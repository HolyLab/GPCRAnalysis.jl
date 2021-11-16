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
        sum(!=(Residue('-')), v)
    end
    mask = nresidues .> minres
    filtersequences!(msa, mask)
end

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
SequenceMapping(seq::AnnotatedAlignedSequence) = SequenceMapping(getsequencemapping(seq))

Base.size(sm::SequenceMapping) = size(sm.seqmap)
Base.getindex(sm::SequenceMapping, i::Int) = sm.seqmap[i]
Base.setindex!(sm::SequenceMapping, v, i::Int) = sm.seqmap[i] = v
Base.similar(::SequenceMapping, ::Type{Int}, dims::Dims{1}) = SequenceMapping(Vector{Int}(undef, dims))

# We restrict the following to `Vector` to prevent ambiguities
Base.getindex(fullseqvec::Vector{T}, sm::SequenceMapping) where T = Union{T,Nothing}[i == 0 ? nothing : fullseqvec[i] for i in sm.seqmap]

function Base.setindex!(fullseqvec::Vector, colvalues::AbstractVector, sm::SequenceMapping)
    idxnz = sm.seqmap .> 0
    fullseqvec[sm.seqmap[idxnz]] = colvalues[idxnz]
end

"""
    chimerax_script(filename, uprot_list, msa::AnnotatedMultipleSequenceAlignment, colidxs;
                    dir=pwd(), chain_transparency=80, styles=Dict{Int,String}(), extras=String[])

Create a [chimerax]() visualization script with name `filename`. `uprot_list` is a list of UniProtX names that you
want to visualize.  The AlphaFold PDB files for these proteins should be `dir`. `msa` is a Multiple Sequence alignment
and `colidxs` specifies the column indices in `msa` that you'd like to visualize.

`chain_transparency` sets the transparency on the ribbon diagrams (0 = not transparent).
`styles` can be used to affect the display, e.g., `Dict(k => "@SD sphere")` would cause methionines at column index `k`
to be displayed with the sulfur in sphere mode. `extras` can be used to hand-specify a number of additional commands;
this can be useful if, for example, the `msa` has occasional misalignments.

# Examples

Suppose you have the `msa` for rhodopsin (mouse: P15409), then:

```julia
chimerax_script("myscript.cxc", ["P15409"], msa, [i1, i2, i3])
```

where `i1` through `i3` are column-indices in the `msa` that you'd like to view.
"""
function chimerax_script(filename, uprot_list, msa::AnnotatedMultipleSequenceAlignment, colidxs;
                         dir=pwd(), chain_transparency=80, styles=Dict{Int,String}(), extras=String[])
    open(filename, "w") do io
        for (i, p) in enumerate(uprot_list)
            afname = alphafoldfile(p)
            println(io, "open ", joinpath(dir, afname))
            seqidx = findfirst(str -> startswith(str, p), sequencenames(msa))
            sm = getsequencemapping(msa, seqidx)
            for c in colidxs
                ridx = sm[c]
                if iszero(ridx)
                    @warn "column $c not set in $p"
                    continue
                end
                println(io, "show #", i, " :", ridx)
                style = get(styles, c, nothing)
                if style !== nothing
                    println(io, "style #", i, " :", ridx, " ", style)
                end
            end
        end
        for ex in extras
            println(io, ex)
        end
        println(io, "transparency #1-", length(uprot_list), ' ', chain_transparency, " target c")
        if !any(str -> startswith(str, "align"), extras)
            println(io, "matchmaker #2-", length(uprot_list), " to #1")
        end
    end
end
