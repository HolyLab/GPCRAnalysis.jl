module GPCRAnalysisBioStockholmExt

using GPCRAnalysis
using BioStockholm
using BioStockholm: OrderedDict   # from OrderedCollections.jl

function conscols(msa::MSA)
    if length(msa.GR) == 1
        # Fast-path: use the reference sequence
        key, _ = only(msa.GR)
        s = msa.seq[key]
        return findall(s) do c
            c == '-' || isuppercase(c)
        end
    end
    # Slow path: check each sequence, find all that have at least one uppercase in that column
    keep = falses(length(msa.GC["seq_cons"]))
    for (_, s) in msa.seq
        keep .|= isuppercase.(s)
    end
    return findall(keep)
end

# Low-level API implementation
GPCRAnalysis.sequenceindexes(msa::MSA, i::Int) = GPCRAnalysis.sequenceindexes(msa::MSA, MSACode(GPCRAnalysis.sequencekeys(msa)[i]))
function GPCRAnalysis.sequenceindexes(msa::MSA, key::MSACode)
    # seq = GPCRAnalysis.msasequence(msa, key)
    seq = msa.seq[String(key)]
    offset = findfirst(!=('.'), seq)
    filled = [r != '-' for r in seq]
    cf = cumsum(filled)
    keepcols = conscols(msa)
    m = match(r"/(\d+)-(\d+)$", String(key))
    if m !== nothing
        start, stop = parse.(Int, m.captures)
        Δ = start - offset
        return (filled .* (cf .+ Δ))[keepcols]
    end
    return (filled .* cf)[keepcols]
end
GPCRAnalysis.sequencekeys(msa::MSA) = collect(keys(msa.seq))
GPCRAnalysis.msasequence(msa::MSA, key::MSACode) = msa.seq[String(key)][conscols(msa)]
GPCRAnalysis.msasequence(msa::MSA, key::AbstractString) = GPCRAnalysis.msasequence(msa, MSACode(key))
function GPCRAnalysis.residuematrix(msa::MSA)
    keepcols = conscols(msa)
    # keepcols = Colon()
    reduce(vcat, [permutedims(seq[keepcols]) for (_, seq) in msa.seq])
end
GPCRAnalysis.subseqs(msa::MSA{T}, rowmask::AbstractVector{Bool}) where T = MSA{T}(OrderedDict(pr for (pr, keep) in zip(msa.seq, rowmask) if keep), msa.GF, OrderedDict(pr for (pr, keep) in zip(msa.GS, rowmask) if keep), msa.GC, msa.GR)
function GPCRAnalysis.subseqs!(msa::MSA, rowmask::AbstractVector{Bool})
    for ((key, _), keep) in zip(msa.seq, rowmask)
        if !keep
            delete!(msa.seq, key)
            delete!(msa.GS, key)
        end
    end
    return msa
end
GPCRAnalysis.columnindexes(msa::BioStockholm.MSA) = conscols(msa)

Base.getindex(msa::MSA, seqname::MSACode) = msa.seq[seqname.name][conscols(msa)]
Base.getindex(msa::MSA, seqname::AccessionCode) = msa[MSACode(msa, seqname)]


function GPCRAnalysis.AccessionCode(msa::MSA, seqname::AbstractString)
    AccessionCode(split(msa.GS[seqname]["AC"], '.')[1])
end
GPCRAnalysis.AccessionCode(msa::MSA, seqname::MSACode) = AccessionCode(msa, seqname.name)
GPCRAnalysis.AccessionCode(::MSA, seqname::AccessionCode) = seqname

function GPCRAnalysis.MSACode(msa::MSA, accession::AbstractString)
    acs = [split(ac["AC"], '.')[1] for (_, ac) in msa.GS]
    i = findfirst(==(accession), acs)
    return MSACode(GPCRAnalysis.sequencekeys(msa)[i])
end
GPCRAnalysis.MSACode(msa::MSA, accession::AccessionCode) = MSACode(msa, accession.name)
GPCRAnalysis.MSACode(::MSA, accession::MSACode) = accession

end
