"""
    getchain(filename::AbstractString; model="1", chain="A")

Read a PDB file `filename` and extract the specified chain.
"""
getchain(filename::AbstractString; model="1", chain="A") = read(filename, PDBFile; model, chain)

function validate_seq_residues(seq::AnnotatedAlignedSequence, chain::AbstractVector{PDBResidue})
    for (i, r) in zip(getsequencemapping(seq), seq)
        (r == GAP || r == XAA) && continue
        res = three2residue(chain[i].id.name)
        res == r || return false
    end
    return true
end

function countitems(list)
    count = Dict{eltype(list),Int}()
    for item in list
        count[item] = get(count, item, 0) + 1
    end
    return count
end

function findall_subseq(subseq, seq)
    starts = Int[]
    n = length(subseq)
    for i in eachindex(seq)
        if seq[i] == subseq[begin] && i + n - 1 <= lastindex(seq)
            issame = true
            for (ss, s) in zip(view(seq, i:i+n-1), subseq)
                issame &= ss == s
            end
            if issame
                push!(starts, i)
            end
        end
    end
    return starts
end

