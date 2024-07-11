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
    getchain(filename::AbstractString; model="1", chain="A")

Read a PDB file `filename` and extract the specified chain.
"""
getchain(filename::AbstractString; model="1", chain="A") = read_file(filename, PDBFile; model, chain)

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

alpha_nitrogen(r) = r.atoms[findfirst(a -> a.atom == "N", r.atoms)]
alpha_carbon(r) = r.atoms[findfirst(a -> a.atom == "CA", r.atoms)]