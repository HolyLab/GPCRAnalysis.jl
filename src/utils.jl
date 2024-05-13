alphafoldfile(uniprotXname; version=4) = "AF-$uniprotXname-F1-model_v$version.pdb"

"""
    fns = alphafoldfiles(dirname=pwd())

Return the latest version of all AlphaFold files in `dirname`.
"""
function alphafoldfiles(dirname=pwd())
    latest = Dict{String,Int}()
    for fn in readdir(dirname)
        m = match(rex_alphafold_pdbs, fn)
        m === nothing && continue
        access_code, v = m.captures
        v = parse(Int, v)
        lv = get(latest, access_code, 0)
        if v > lv
            latest[access_code] = v
        end
    end
    return sort!([alphafoldfile(access_code; version=v) for (access_code, v) in latest])
end

"""
    msacode2structfile = alphafoldfiles(msa::AnnotatedMultipleSequenceAlignment, dirname=pwd())

Return a dictionary mapping `MSACode`s to the corresponding AlphaFold structure files.
"""
function alphafoldfiles(msa::AnnotatedMultipleSequenceAlignment, dirname=pwd())
    afs = alphafoldfiles(dirname)
    accesscode2idx = Dict{AccessionCode,Int}()
    for (i, af) in pairs(afs)
        ac = AccessionCode(match(rex_alphafold_pdbs, af).captures[1])
        accesscode2idx[ac] = i
    end
    msacode2structfile = Dict{MSACode,String}()
    for name in sequencenames(msa)
        ac = AccessionCode(msa, name)
        if haskey(accesscode2idx, ac)
            msacode2structfile[MSACode(name)] = afs[accesscode2idx[ac]]
        end
    end
    return msacode2structfile
end

function query_alphafold_latest(accession_code)
    resp = HTTP.get("https://alphafold.com/api/prediction/$accession_code")
    if resp.status == 200
        j = JSON.parse(String(resp.body))[1]
        return j["pdbUrl"]
    end
    return nothing
end


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
function download_alphafolds(msa::AbstractMultipleSequenceAlignment; dirname=pwd(), maxversion=nothing, kwargs...)
    maxversion === nothing || @warn "`download_alphafolds`: `maxversion` kwarg has no effect and is deprecated" maxlog=1
    @showprogress 1 "Downloading AlphaFold files..." for name in sequencenames(msa)
        uname = accession_code(msa, name)
        url = query_alphafold_latest(uname)
        url === nothing && continue
        fn = split(url, '/')[end]
        path = joinpath(dirname, fn)
        if !isfile(path)
            Downloads.download(url, path)
        end
        if !validate_seq_residues(getsequence(msa, name), getchain(path))
            @warn "Residues in $path do not match those in the sequence $name, removing PDB file"
            rm(path)
        end
    end
end

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

function align_to_axes(strct)
    coords = coordinatesmatrix(strct)
    coords = coords .- mean(coords, dims=1)
    tform = GaussianMixtureAlignment.inertial_transforms(coords')[1]
    tform = AffineMap(Matrix(tform.linear),[tform.translation...])
    return change_coordinates(strct, tform(coords')')
end
