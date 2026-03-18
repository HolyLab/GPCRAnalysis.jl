formats = Dict(:json => "application/json",
               :fasta => "text/x-fasta",
               :gff => "application/gff",
               :gtf => "application/gtf",
               :gb => "application/genbank",
               :embl => "application/embl",
               :xml => "application/xml",
               :txt => "text/plain")

"""
    result = query_ebi_proteins(id; format=:json)

Query the EBI Proteins API for a protein with the specified `id`, which must be the Uniprot accession code.
You can also supply several proteins as a comma-separated list.

`result` is a JSON3 object with many fields.
"""
function query_ebi_proteins(id::AbstractString; size=count(==(','), id)+1, format::Symbol=:json)
    resp = HTTP.get("https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=$size&accession=$id", ["Accept" => formats[format]]; status_exception=false)
    if resp.status == 200
        return mktemp() do tmppath, tmpio
            body = String(resp.body)
            write(tmpio, body)
            flush(tmpio)
            uncompr_body = GZip.open(tmppath) do iogz
                read(iogz, String)
            end
            return format===:json ? JSON3.read(uncompr_body) : uncompr_body
        end
    end
    return nothing
end
query_ebi_proteins(id::AbstractVector{<:AbstractString}; kwargs...) =
    query_ebi_proteins(join(id, ','); kwargs...)

"""
    ntms = query_tm_count(id)
    ntms = query_tm_count(ids)

Query the EBI Proteins API for the number of annotated transmembrane (TRANSMEM) segments.

The single-accession form returns an `Int` (or `nothing` if the entry is not found). The
vector form returns a `Dict{String,Int}` mapping each returned accession code to its count.
Queries are batched at 50 IDs per request.

Useful for choosing an `id_for_tms` for [`align_family`](@ref): a good reference receptor
should have its TM segments fully annotated (e.g., 7 for a canonical GPCR).
"""
function query_tm_count(ids::AbstractVector{<:AbstractString}; batchsize=50)
    result = Dict{String,Int}()
    for i in firstindex(ids):batchsize:lastindex(ids)
        batch = ids[i:min(i+batchsize-1, lastindex(ids))]
        response = query_ebi_proteins(batch)
        response === nothing && continue
        for r in response
            result[String(r["accession"])] = count(f -> f["type"] == "TRANSMEM", r["features"])
        end
    end
    return result
end

function query_tm_count(id::AbstractString)
    result = query_tm_count([id])
    return get(result, id, nothing)
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
