# Support output of TM-align: https://zhanggroup.org/TM-align/
# Zhang, Yang, and Jeffrey Skolnick. "TM-align: a protein structure alignment algorithm based on the TM-score."
# Nucleic acids research 33.7 (2005): 2302-2309.

struct MapAlign
    s2a::Vector{Int}    # sequence index => alignment index
    a2s::Vector{Int}    # alignment index => sequence index
end

struct StructAlign
    m1::MapAlign
    m2::MapAlign

    function StructAlign(m1::MapAlign, m2::MapAlign)
        axes(m1.a2s) == axes(m2.a2s) || error("the two maps disagree")
        return new(m1, m2)
    end
end

function parse_tma(filename)
    lines = readlines(filename)
    filter!(!isempty, lines)
    idx = findfirst(line -> occursin("denotes", line), lines)  # check for the header
    if idx !== nothing
        lines = lines[idx+1:end]
    end
    length(lines) == 3 || error("expecting 3 lines, got $(length(lines))")
    seq1, seq2 = collect(lines[1]), collect(lines[3])
    quality = collect(lines[2])
    append!(quality, fill('-', length(seq1) - length(quality)))
    return seq1, seq2, quality
end

const alignsentinel = -1

function MapAlign(s::ChainLike, a::AbstractVector, quality)
    s2a, a2s = fill(alignsentinel, length(s)), Int[]  # we only retain the high-quality matches
    sj, state = iterate(s)
    j = 0
    for (r, q) in zip(a, quality)
        isgap(r) && continue
        j += 1
        Char(r) == three2char(resname(sj)) || error("at position $j, residue was $(resname(sj)) but alignment was $r")
        if q == ':'
            push!(a2s, j)
            s2a[j] = length(a2s)
        end
        ret = iterate(s, state)
        sj, state = ret === nothing ? (nothing, nothing) : ret   # if we've reached s[end], continuing should be an error
    end
    return MapAlign(s2a, a2s)
end

"""
    StructAlign(struct1::ChainLike, struct2::ChainLike, filename::AbstractString)

Create a structure-based alignment between `struct1` and `struct2`. `filename` is the name of the TM-align "results" file
(e.g., https://zhanggroup.org//TM-align/example/873772.html).

See also [`residueindex`](@ref).
"""
StructAlign(struct1::ChainLike, struct2::ChainLike, filename::AbstractString) =
    StructAlign(struct1, struct2, parse_tma(filename)...)

function StructAlign(struct1::ChainLike, struct2::ChainLike,
                     align1::AbstractVector{Char}, align2::AbstractVector{Char},
                     quality)
    StructAlign(MapAlign(struct1, align1, quality), MapAlign(struct2, align2, quality))
end


"""
    residueindex(sa::StructAlign, nothing, idx2)
    residueindex(sa::StructAlign, nothing, idx2, ±1)

Calculate the residue index in structure 1 corresponding to residue `idx2` in structure 2.
The second form allows you to find the nearest corresponding residue in the forward (+1) or reverse (-1)
directions, if `idx2` is not among the mapped residues.
"""
residueindex(sa::StructAlign, ::Nothing, idx2::Integer) = sa.m1.a2s[sa.m2.s2a[idx2]]
function residueindex(sa::StructAlign, ::Nothing, idx2::Integer, dir)
    dir ∈ (-1, 1) || error("dir must be ±1")
    i2 = sa.m2.s2a[idx2]
    while i2 == alignsentinel
        idx2 += dir
        i2 = sa.m2.s2a[idx2]
    end
     return sa.m1.a2s[i2]
end
ismapped(sa::StructAlign, ::Nothing, idx2::Integer) = sa.m2.s2a[idx2] != alignsentinel

"""
    residueindex(sa::StructAlign, idx1, nothing)
    residueindex(sa::StructAlign, idx1, nothing, ±1)

Calculate the residue index in structure 2 corresponding to residue `idx1` in structure 1.
The second form allows you to find the nearest corresponding residue in the forward (+1) or reverse (-1)
directions, if `idx1` is not among the mapped residues.
"""
residueindex(sa::StructAlign, idx1::Integer, ::Nothing) = sa.m2.a2s[sa.m1.s2a[idx1]]
function residueindex(sa::StructAlign, idx1::Integer, ::Nothing, dir)
    dir ∈ (-1, 1) || error("dir must be ±1")
    i1 = sa.m1.s2a[idx1]
    while i1 == alignsentinel
        idx1 += dir
        i1 = sa.m1.s2a[idx1]
    end
     return sa.m2.a2s[i1]
end
ismapped(sa::StructAlign, idx1::Integer, ::Nothing) = sa.m1.s2a[idx1] != alignsentinel
