alphafoldfile(uniprotXname) = "AF-$uniprotXname-F1-model_v1.pdb"

"""
    try_download_alphafold(uniprotXname, path=alphafoldfile(uniprotXname))

Attempt to download an [AlphaFold](https://alphafold.com/) structure.
Returns `nothing` if no entry corresponding to `uniprotXname` exists; otherwise
it return `path`, the pathname of the saved file.
"""
function try_download_alphafold(uniprotXname::AbstractString, path::AbstractString=alphafoldfile(uniprotXname))
    fn = alphafoldfile(uniprotXname)
    try
        Downloads.download("https://alphafold.ebi.ac.uk/files/$fn", path)
    catch err
        if isa(err, RequestError) && err.response.status == 404
            return nothing
        else
            throw(err)
        end
    end
    return path
end

"""
    download_alphafolds(msa::AbstractMultipleSequenceAlignment; dirname=pwd())

Download all available [AlphaFold](https://alphafold.com/) structures for the sequences in `msa`.
Missing entries are silently skipped.
"""
function download_alphafolds(msa::AbstractMultipleSequenceAlignment; dirname=pwd())
    for name in sequencenames(msa)
        uname = uniprotX(name)
        path = joinpath(dirname, alphafoldfile(uname))
        if !isfile(path)
            try_download_alphafold(uname, path)
        end
    end
end

"""
    getchain(filename::AbstractString; model="1", chain="A")

Read a PDB file `filename` and extract the specified chain.
"""
getchain(filename::AbstractString; model="1", chain="A") = read(filename, PDBFile; model, chain)

function countitems(list)
    count = Dict{eltype(list),Int}()
    for item in list
        count[item] = get(count, item, 0) + 1
    end
    return count
end
