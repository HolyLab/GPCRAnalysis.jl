alphafoldfile(uniprotXname; version=4) = "AF-$uniprotXname-F1-model_v$version.pdb"

"""
    try_download_alphafold(uniprotXname, path=alphafoldfile(uniprotXname))

Attempt to download an [AlphaFold](https://alphafold.com/) structure.
Returns `nothing` if no entry corresponding to `uniprotXname` exists; otherwise
it return `path`, the pathname of the saved file.
"""
function try_download_alphafold(uniprotXname::AbstractString, path::AbstractString=alphafoldfile(uniprotXname); kwargs...)
    fn = alphafoldfile(uniprotXname; kwargs...)
    isfile(path) && return path
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
function download_alphafolds(msa::AbstractMultipleSequenceAlignment; dirname=pwd(), maxversion=4)
    @showprogress 1 "Downloading AlphaFold files..." for name in sequencenames(msa)
        uname = uniprotX(name)
        path = joinpath(dirname, alphafoldfile(uname))
        if !isfile(path)
            for version = maxversion:-1:1
                fn = try_download_alphafold(uname, path; version)
                fn !== nothing && break
            end
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
