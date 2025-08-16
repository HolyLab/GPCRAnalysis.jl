# Generates 2 captures, one for the uniprotXname and the other for the version
const rex_alphafold_pdbs = r"AF-([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})-F1-model_v(\d+).(?:pdb|cif|bcif)"
# Make a regex for a specific uniprotXname (single capture for the version)
regex_alphafold_pdb(uniprotXname) = Regex("AF-$uniprotXname-F1-model_v(\\d+).(?:pdb|cif|bcif)")

alphafoldfilename(uniprotXname; ext="pdb", version=4) = "AF-$uniprotXname-F1-model_v$version.$ext"

"""
    fns = alphafoldfile(uniprotXname, dirname=pwd(); join=false)

Return the latest version of the AlphaFold file for `uniprotXname` in `dirname`.
If `join` is `true`, then the full path is returned.
"""
function alphafoldfile(uniprotXname::AbstractString, dirname=pwd(); join::Bool=false)
    rex = regex_alphafold_pdb(uniprotXname)
    fnv, lv = "", 0
    for fn in readdir(dirname)
        m = match(rex, fn)
        m === nothing && continue
        v = parse(Int, only(m.captures))
        if v > lv
            fnv, lv = fn, v
        end
    end
    lv == 0 && return nothing
    return join ? joinpath(dirname, fnv) : fnv
end

"""
    fns = alphafoldfiles(dirname=pwd(); join=false)

Return the latest version of all AlphaFold files in `dirname`.
If `join` is `true`, then the full paths are returned.
"""
function alphafoldfiles(dirname::AbstractString=pwd(); join::Bool=false)
    latest = Dict{String,Int}()
    latestfn = Dict{String,String}()
    for fn in readdir(dirname)
        m = match(rex_alphafold_pdbs, fn)
        m === nothing && continue
        access_code, v = m.captures
        v = parse(Int, v)
        lv = get(latest, access_code, 0)
        if v > lv
            latest[access_code] = v
            latestfn[access_code] = fn
        end
    end
    fns = sort!(collect(values(latestfn)))
    return join ? [joinpath(dirname, fn) for fn in fns] : fns
end

"""
    msacode2structfile = alphafoldfiles(msa, dirname=pwd())

Return a dictionary mapping `MSACode`s to the corresponding AlphaFold structure files.
"""
function alphafoldfiles(msa, dirname::AbstractString=pwd(); join::Bool=false)
    afs = alphafoldfiles(dirname)
    accesscode2idx = Dict{AccessionCode,Int}()
    for (i, af) in pairs(afs)
        ac = AccessionCode(match(rex_alphafold_pdbs, af).captures[1])
        accesscode2idx[ac] = i
    end
    msacode2structfile = Dict{MSACode,String}()
    for name in sequencekeys(msa)
        ac = AccessionCode(msa, name)
        if haskey(accesscode2idx, ac)
            fn = afs[accesscode2idx[ac]]
            msacode2structfile[MSACode(name)] = join ? joinpath(dirname, fn) : fn
        end
    end
    return msacode2structfile
end

"""
    url = query_alphafold_latest(uniprotXname; format="cif")

Query the [AlphaFold](https://alphafold.com/) API for the latest structure of `uniprotXname`.
`format` should be "cif", "pdb", or "bcif".
"""
function query_alphafold_latest(uniprotXname::AbstractString; format="cif")
    resp = HTTP.get("https://alphafold.com/api/prediction/$uniprotXname?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94", ["Accept" => "application/json"]; status_exception = false)
    if resp.status == 200
        j = JSON3.read(String(resp.body))[1]
        return j["$(format)Url"]
    end
    return nothing
end
query_alphafold_latest(uniprotXname; kwargs...) = query_alphafold_latest(String(uniprotXname)::String; kwargs...)

"""
    try_download_alphafold(uniprotXname, path=alphafoldfilename(uniprotXname); version=4)

Attempt to download an [AlphaFold](https://alphafold.com/) structure.
Returns `nothing` if no entry corresponding to `uniprotXname` exists; otherwise
it returns `path`, the pathname of the saved file.

In general, a better approach is to use [`download_alphafolds`](@ref) for multiple proteins,
or [`query_alphafold_latest`](@ref) combined with `Downloads.download` for a single protein.
"""
function try_download_alphafold(uniprotXname::AbstractString, path::AbstractString=alphafoldfilename(uniprotXname); kwargs...)
    ext = splitext(path)[2][2:end]  # remove the leading dot
    fn = alphafoldfilename(uniprotXname; ext, kwargs...)
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
    download_alphafolds(msa; dirname=pwd())
    download_alphafolds(ids; dirname=pwd())

Download all available [AlphaFold](https://alphafold.com/) structures for the
sequences in `msa`. Missing entries are silently skipped.

If an `msa` is provided, each downloaded PDB file is checked to ensure that the
residues in the MSA sequence match those in the PDB file. If they do not match,
the PDB file is removed.
"""
function download_alphafolds(msa; dirname=pwd(), maxversion=nothing, kwargs...)
    maxversion === nothing || @warn "`download_alphafolds`: `maxversion` kwarg has no effect and is deprecated" maxlog=1
    @showprogress 1 "Downloading AlphaFold files..." for name in sequencekeys(msa)
        uname = AccessionCode(msa, name)
        url = query_alphafold_latest(uname; kwargs...)
        url === nothing && continue
        fn = split(url, '/')[end]
        path = joinpath(dirname, fn)
        if !isfile(path)
            Downloads.download(url, path)
        end
        if !validate_seq_residues(msasequence(msa, name), getchain(path))
            @warn "Residues in $path do not match those in the sequence $name, removing structure file"
            rm(path)
        end
    end
end

function download_alphafolds(ids::AbstractVector{<:AbstractString}; dirname=pwd(), kwargs...)
    @showprogress 1 "Downloading AlphaFold files..." for uname in ids
        url = query_alphafold_latest(uname; kwargs...)
        url === nothing && continue
        fn = split(url, '/')[end]
        path = joinpath(dirname, fn)
        if !isfile(path)
            Downloads.download(url, path)
        end
    end
end

const pLDDTcolors = [  # see https://github.com/sokrypton/ColabFold?tab=readme-ov-file#faq
    RGB{N0f8}(0.051, 0.341, 0.827),
    RGB{N0f8}(0.416, 0.796, 0.945),
    RGB{N0f8}(0.996, 0.851, 0.212),
    RGB{N0f8}(0.992, 0.490, 0.302),
]
function pLDDTcolor(pLDDT::Real)
    if 0 <= pLDDT < 50
        return pLDDTcolors[4]
    elseif pLDDT < 70
        return pLDDTcolors[3]
    elseif pLDDT < 90
        return pLDDTcolors[2]
    elseif pLDDT <= 100
        return pLDDTcolors[1]
    else
        throw(ArgumentError("pLDDT score must be between 0 and 100, got $pLDDT"))
    end
end

"""
    pLDDTcolor(r::Residue)
    pLDDTcolor(score::Real)

Return the color corresponding to the pLDDT score (a measure of confidence) of a residue.
"""
pLDDTcolor(r::Residue) = pLDDTcolor(pLDDT(r))
pLDDTcolor(a::Atom) = pLDDTcolor(pLDDT(a))

pLDDT(r::Residue) = only(unique(pLDDT, r))
pLDDT(a::Atom) = tempfactor(a)
