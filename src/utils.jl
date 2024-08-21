"""
    result = query_ebi_proteins(id)

Query the EBI Proteins API for a protein with the specified `id`, which must be the Uniprot accession code.
You can also supply several proteins as a comma-separated list.

`result` is a JSON3 object with many fields.
"""
function query_ebi_proteins(id; size=count(==(','), id)+1)
    resp = HTTP.get("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=$size&accession=$id", ["Accept" => "application/json"])
    if resp.status == 200
        return mktemp() do tmppath, tmpio
            body = String(resp.body)
            write(tmpio, body)
            flush(tmpio)
            uncompr_body = GZip.open(tmppath) do iogz
                read(iogz, String)
            end
            return JSON3.read(uncompr_body)
        end
    end
    return nothing
end

"""
    result = query_ncbi(id)

Query the NCBI API for a gene with the specified `id`, which must be the gene accession code.
`result` is a JSON3 object with many fields.
"""
function query_ncbi(id)
    resp = HTTP.get("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/accession/$id?table_fields=gene-id&table_fields=gene-type&table_fields=description",
                    ["Accept" => "application/json"])
    if resp.status == 200
        return JSON3.read(String(resp.body))
    end
    return nothing
end

"""
    getchain(filename::AbstractString; model=1, chain="A")

Read a PDB file `filename` and extract the specified chain.
"""
getchain(filename::AbstractString; model=1, chain="A") = read(filename, PDBFormat)[model][chain]

function validate_seq_residues(seq::AnnotatedAlignedSequence, chain)
    for (i, r) in zip(getsequencemapping(seq), seq)
        (r == GAP || r == XAA) && continue
        res = three2residue(String(resname(chain[i])))
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

MIToS.PDB.ishydrophobic(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in MIToS.PDB._hydrophobic
MIToS.PDB.isaromatic(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in MIToS.PDB._aromatic
MIToS.PDB.iscationic(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in MIToS.PDB._cationic
MIToS.PDB.isanionic(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in MIToS.PDB._anionic
MIToS.PDB.ishbonddonor(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in keys(MIToS.PDB._hbond_donor)
MIToS.PDB.ishbondacceptor(a::AbstractAtom, rname::AbstractString) = (rname, atomname(a)) in keys(MIToS.PDB._hbond_acceptor)
