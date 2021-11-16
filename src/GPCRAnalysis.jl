module GPCRAnalysis

using Downloads
using Statistics

using MIToS
using MIToS.Pfam
using MIToS.MSA
using MIToS.PDB

using MultivariateStats
using Distances

export @res_str

export SequenceMapping
export species, uniprotX
export try_download_alphafold, download_alphafolds, getchain, findall_subseq
export filter_species!, filter_long!, chimerax_script
export StructAlign, residueindex
export project_sequences, columnwise_entropy, align, residue_centroid, residue_centroid_matrix, mapclosest, chargelocations, positive_locations, negative_locations
export BWScheme, lookupbw

include("naming_conventions.jl")
include("utils.jl")
include("msa.jl")
include("analyze.jl")
include("tmalign.jl")
include("ballesteros_weinstein.jl")

end
