three2char(resname::AbstractString) = Char(threeletter_to_aa[String(resname)])
three2char(r::AbstractResidue) = three2char(resname(r))

"""
    getchain(filename::AbstractString; model=1, chain="A")

Read a PDB or mmCIF file `filename` and extract the specified chain.
"""
getchain(filename::AbstractString; model=1, chain="A") =
    endswith(filename, ".pdb") ? read(filename, PDBFormat)[model][chain] :
    endswith(filename, ".cif") ? read(filename, MMCIFFormat)[model][chain] :
    throw(ArgumentError("Unknown format for $filename"))

"""
    writechain(filename::AbstractString, chain::ChainLike)

Write the specified `chain` to a PDB or mmCIF file `filename`.
"""
writechain(filename::AbstractString, chain::ChainLike) =
    endswith(filename, ".pdb") ? writepdb(filename, chain) :
    endswith(filename, ".cif") ? writemmcif(filename, chain) :
    throw(ArgumentError("Unknown format for $filename"))

"""
    validate_seq_residues(msaseq, chain)

Return `true` if the residues in `msaseq` match those in `chain`, ignoring gaps and unknown residues.
"""
function validate_seq_residues(msaseq, chain)
    try
        for (i, r) in zip(sequenceindexes(msaseq), msaseq)
            (isgap(r) || isunknown(r)) && continue
            res = three2char(String(resname(chain[i])))
            res == Char(r) || return false
        end
    catch
        @warn "Error validating sequence residues"
        return false
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

alpha_nitrogen(r) = firstmatch(a -> atomname(a) == "N", r)
alpha_carbon(r) = firstmatch(a -> atomname(a) == "CA", r)

function firstmatch(f, iter)
    for x in iter
        if f(x)
            return x
        end
    end
    return nothing
end

function skipnothing(v::AbstractVector{Union{T,Nothing}}) where T
    return T[x for x in v if x !== nothing]
end
