module GPCRAnalysisMIToSExt

using GPCRAnalysis
using Downloads
using BioStructures
using ProgressMeter

using GPCRAnalysis: ChainLike, ResidueLike, StructureLike, _entropy, validate_seq_residues, rex_alphafold_pdbs

using MIToS: MIToS, Pfam, MSA
using MIToS.MSA: AbstractMultipleSequenceAlignment, AnnotatedAlignedSequence, AnnotatedMultipleSequenceAlignment,
                 ReducedAlphabet, ResidueAlphabet, GAP, XAA
using MIToS.MSA: getsequence, getannotsequence, getsequencemapping, getresidues, three2residue, sequencenames,
                 filtersequences, filtersequences!, percentsimilarity, getcolumnmapping


# Low-level API implementation
GPCRAnalysis.sequenceindexes(msaseq::AnnotatedAlignedSequence) = getsequencemapping(msaseq)
GPCRAnalysis.sequenceindexes(msa::AbstractMultipleSequenceAlignment, i::Int) = getsequencemapping(msa, i)
GPCRAnalysis.sequenceindexes(msa::AbstractMultipleSequenceAlignment, key::AbstractString) = getsequencemapping(msa, key)
GPCRAnalysis.sequenceindexes(msa::AbstractMultipleSequenceAlignment, key::MSACode) = sequenceindexes(msa, String(key))
GPCRAnalysis.isgap(res::MSA.Residue) = res == GAP
GPCRAnalysis.isunknown(res::MSA.Residue) = res == XAA
GPCRAnalysis.sequencekeys(msa::AbstractMultipleSequenceAlignment) = sequencenames(msa)
GPCRAnalysis.msasequence(msa::AbstractMultipleSequenceAlignment, key::AbstractString) = getsequence(msa, key)
GPCRAnalysis.msasequence(msa::AbstractMultipleSequenceAlignment, key::MSACode) = msasequence(msa, String(key))
GPCRAnalysis.residuematrix(msa::AbstractMultipleSequenceAlignment) = getresidues(msa)
GPCRAnalysis.subseqs(msa::AbstractMultipleSequenceAlignment, rowmask)  = filtersequences(msa, rowmask)
GPCRAnalysis.subseqs!(msa::AbstractMultipleSequenceAlignment, rowmask) = filtersequences!(msa, rowmask)
GPCRAnalysis.percent_similarity(msa::AbstractMultipleSequenceAlignment) = percentsimilarity(msa)
GPCRAnalysis.columnindexes(msa::MSA.AbstractMultipleSequenceAlignment) = getcolumnmapping(msa)

Base.getindex(msa::AbstractMultipleSequenceAlignment, seqname::MSACode) = getsequence(msa, seqname.name)
Base.getindex(msa::AbstractMultipleSequenceAlignment, seqname::AccessionCode) = getsequence(msa, MSACode(msa, seqname).name)

function GPCRAnalysis.AccessionCode(msa::AnnotatedMultipleSequenceAlignment, seqname::AbstractString)
    AccessionCode(uniprotX(getannotsequence(msa, seqname, "AC", seqname)))
end
GPCRAnalysis.AccessionCode(msa::AnnotatedMultipleSequenceAlignment, seqname::MSACode) = AccessionCode(msa, seqname.name)
GPCRAnalysis.AccessionCode(::AnnotatedMultipleSequenceAlignment, seqname::AccessionCode) = seqname

function GPCRAnalysis.MSACode(msa::AnnotatedMultipleSequenceAlignment, accession::AbstractString)
    seqnames = sequencenames(msa)
    return MSACode(seqnames[findfirst(x -> AccessionCode(msa, x).name == accession, seqnames)])
end
GPCRAnalysis.MSACode(msa::AnnotatedMultipleSequenceAlignment, accession::AccessionCode) = MSACode(msa, accession.name)
GPCRAnalysis.MSACode(::AnnotatedMultipleSequenceAlignment, accession::MSACode) = accession

GPCRAnalysis.SequenceMapping(seq::AnnotatedAlignedSequence) = SequenceMapping(getsequencemapping(seq))

# Move this to MIToS?
if !hasmethod(getsequencemapping, Tuple{AnnotatedAlignedSequence})
    function MIToS.MSA.getsequencemapping(seq::AnnotatedAlignedSequence)
        getsequencemapping(seq, sequencenames(seq)[1])
    end
    function MIToS.MSA.getsequencemapping(msa::Union{AnnotatedAlignedSequence,AnnotatedMultipleSequenceAlignment}, seq_id::String)
        MIToS.MSA._str2int_mapping(getannotsequence(msa, seq_id, "SeqMap"))
    end
    function MIToS.MSA.getsequencemapping(msa::AnnotatedMultipleSequenceAlignment, seqid::Regex)
        id = findfirst(str -> occursin(seqid, str), sequencenames(msa))
        getsequencemapping(msa, id)
    end
end

const reduced_code = ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")

"""
    columnwise_entropy(msa, aacode = reduced_code)

Call `columnwise_entropy` after mapping each residue through `aacode`.

The default code is `ReducedAlphabet("(AILMV)(NQST)(RHK)(DE)(FWY)CGP")`, which
groups residues into categories hydrophobic, polar, charged, aromatic, and
"special."
"""
GPCRAnalysis.columnwise_entropy(msa::AbstractMultipleSequenceAlignment, aacode::ResidueAlphabet=reduced_code) =
    GPCRAnalysis.columnwise_entropy(r -> aacode[r], msa)

end
