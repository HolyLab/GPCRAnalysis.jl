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
