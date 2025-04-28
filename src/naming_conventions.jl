abstract type NamingConvention end

struct MSACode <: NamingConvention
    name::String
end
struct AccessionCode <: NamingConvention
    name::String
end

Base.String(x::NamingConvention) = x.name
Base.convert(::Type{NC}, x::AbstractString) where NC<:NamingConvention = NC(x)  # but don't support the other direction
Base.:(==)(a::T, b::T) where T<:NamingConvention = a.name == b.name

const rex_resrange = r"/(\d+)-(\d+)$"
function residue_range(name::AbstractString)
    m = match(rex_resrange, name)
    m === nothing && error("no residue range found")
    return parse(Int, m.captures[1]):parse(Int, m.captures[2])
end
function strip_residue_range(name::AbstractString)
    m = match(rex_resrange, name)
    m === nothing && error("no residue range found")
    return name[begin:prevind(name, m.offset)]
end

const rex_uniprot_species = r"_([A-Z0-9]{1,5})($|/)"

"""
    species(name)

Extract the species identifier from a [UniProt "X_Y" entry](https://www.uniprot.org/help/entry_name)
or elaborated variant (e.g., PFAM sequence name).

See also [`uniprotX`](@ref).

# Examples

```jldoctest
julia> species("Q8VGW6_MOUSE/31-308")
"MOUSE"
```
"""
function species(name::AbstractString)
    m = match(rex_uniprot_species, name)
    m === nothing && error("no species found")
    return m.captures[1]::AbstractString
end

const rex_uniprotX_Swiss = r"^([A-Z0-9]{1,5})(?:$|_)"   # Swiss
const rex_uniprot_accession = r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:$|_|\||\.\d+)"

"""
    uniprotX(name)

Extract the [UniProt "X" entry](https://www.uniprot.org/help/entry_name) from an X_Y entry or elaborated variant (e.g., PFAM sequence name).

See also [`species`](@ref).

# Examples

```jldoctest
julia> uniprotX("Q8VGW6_MOUSE/31-308")
"Q8VGW6"
```
"""
function uniprotX(name::AbstractString)::AbstractString
    if startswith(name, "tr|") || startswith(name, "sp|")
        name = name[4:end]
    end
    m = match(rex_uniprotX_Swiss, name)
    m !== nothing && return m.captures[1]
    m = match(rex_uniprot_accession, name)
    m !== nothing && return m.captures[1]
    error("no UniProt X identifier found")
end

"""
    accession_code = query_uniprot_accession(id)

Perform a Uniprot search for `id`, returning the canonical accession code.

# Examples

```
julia> query_uniprot_accession("T2R38_MOUSE")
"Q7TQA6"
```
"""
function query_uniprot_accession(id)
    resp = HTTP.get("https://rest.uniprot.org/uniprotkb/search?query=id:$id&format=json")
    if resp.status == 200
        return mktemp() do tmppath, tmpio
            body = String(resp.body)
            write(tmpio, body)
            flush(tmpio)
            uncompr_body = GZip.open(tmppath) do iogz
                read(iogz, String)
            end
            j = JSON3.read(uncompr_body)
            return j["results"][1]["primaryAccession"]
        end
    end
    return nothing
end

"""
    jobID = GPCRAnalysis.map_uniprot_submit(ids, from="UniProtKB_AC-ID", to="UniProtKB")

Submit a list of `ids` to the Uniprot ID mapping service, to convert from ID convention `from` to `to`.
The jobID can be used to check the status (`map_uniprot_status`) and retrieve the results (`map_uniprot_retrieve`).

# Examples
```jldoctest`
julia> jobID = GPCRAnalysis.map_uniprot_submit(["ENSMUSG00000067064", "ENSMUSG00000057464"], "Ensembl");
```
"""
function map_uniprot_submit(ids::AbstractString, from::AbstractString="UniProtKB_AC-ID", to::AbstractString="UniProtKB")
    resp = HTTP.post("https://rest.uniprot.org/idmapping/run", [],
                     Dict("from" => from, "to" => to, "ids" => ids))
    if resp.status == 200
        return JSON3.read(String(resp.body))["jobId"]
    end
    return nothing
end
map_uniprot_submit(ids::AbstractVector, args...) = map_uniprot_submit(join(ids, ','), args...)

"""
    status = GPCRAnalysis.map_uniprot_status(jobID)

Check the status of a Uniprot ID mapping job. Returns `true` if the results are
ready. Otherwise, returns the status object.
"""
function map_uniprot_status(jobID)
    resp = HTTP.get("https://rest.uniprot.org/idmapping/status/$jobID", ["Accept" => "application/json"]; decompress = true)
    if resp.status == 200
        status = JSON3.read(String(HTTP.decode(resp)))
        haskey(status, "results") && return true
        return status
    end
    return nothing
end

"""
    result = GPCRAnalysis.map_uniprot_retrieve(jobID)

Retrieve the results of a Uniprot ID mapping job.
"""
function map_uniprot_retrieve(jobID)
    resp = HTTP.get("https://rest.uniprot.org/idmapping/stream/$jobID", ["Accept" => "application/json"]; decompress = true)
    if resp.status == 200
        return JSON3.read(String(HTTP.decode(resp)))
    end
    return nothing
end
