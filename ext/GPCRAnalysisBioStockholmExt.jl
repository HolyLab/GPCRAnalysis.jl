module GPCRAnalysisBioStockholmExt

using GPCRAnalysis
using BioStockholm
using BioStockholm: OrderedDict   # from OrderedCollections.jl

function conscols(msa::MSA)
    ss = get(msa.GC, "SS_cons", nothing)
    if ss === nothing
        ss = msa.GC["seq_cons"]
    end
    return findfirst(!=( '.'), ss):findlast(!=( '.'), ss)
end

# Low-level API implementation
# GPCRAnalysis.sequenceindexes(msaseq::AnnotatedAlignedSequence) = getsequencemapping(msaseq)
GPCRAnalysis.sequenceindexes(msa::MSA, i::Int) = GPCRAnalysis.sequenceindexes(msa::MSA, MSACode(GPCRAnalysis.sequencekeys(msa)[i]))
function GPCRAnalysis.sequenceindexes(msa::MSA, key::MSACode)
    seq = GPCRAnalysis.msasequence(msa, key)
    filled = [r ∉ ('-', '.') for r in seq]
    start, stop = parse.(Int, match(r"/(\d+)-(\d+)$", String(key)).captures)
    rng = start:stop
    cf = cumsum(filled)
    # cf[end] == length(rng) || return nothing
    return [filled[i] ? cf[i] + start - 1 : 0 for i in eachindex(filled)]
end
GPCRAnalysis.sequencekeys(msa::MSA) = collect(keys(msa.seq))
GPCRAnalysis.msasequence(msa::MSA, key::MSACode) = msa.seq[String(key)][conscols(msa)]
GPCRAnalysis.msasequence(msa::MSA, key::AbstractString) = GPCRAnalysis.msasequence(msa, MSACode(key))
function GPCRAnalysis.residuematrix(msa::MSA)
    keepcols = conscols(msa)
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


function reduced_alphabet(r::Char)
    if r == '-'
        return 0
    elseif r in ('A','I','L','M','V')
        return 1  # hydrophobic
    elseif r in ('N','Q','S','T')
        return 2  # polar
    elseif r in ('R','H','K')
        return 3  # charged
    elseif r in ('D','E')
        return 4  # charged
    elseif r in ('F','W','Y')
        return 5  # aromatic
    end
    offset = findfirst(==(r), ('C','G','P'))
    offset === nothing && throw(ArgumentError("Unknown residue '$r'"))
    return 5 + offset  # special or unknown
end

GPCRAnalysis.columnwise_entropy(msa) = columnwise_entropy(reduced_alphabet, msa)

function GPCRAnalysis.percent_similarity(f, msa::MSA)
    # This mimics MIToS's implementation
    function pctsim(v1, v2)
        same = l = 0
        for (a, b) in zip(v1, v2)
            a == b == 0 && continue  # skip gaps
            same += a == b
            l += 1
        end
        return 100 * same / l
    end

    M = f.(GPCRAnalysis.residuematrix(msa))
    n = size(M, 1)
    S = zeros(Float64, n, n)
    for i in 1:n
        for j in i:n
            S[i, j] = pctsim(M[i, :], M[j, :])
            S[j, i] = S[i, j]
        end
    end
    return S
end
GPCRAnalysis.percent_similarity(msa::MSA) = GPCRAnalysis.percent_similarity(reduced_alphabet, msa)


end
