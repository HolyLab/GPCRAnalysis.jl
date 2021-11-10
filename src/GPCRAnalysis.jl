module GPCRAnalysis

using Downloads
using Statistics

using MIToS
using MIToS.Pfam
using MIToS.MSA
using MIToS.PDB

using MultivariateStats

export @res_str

export SequenceMapping
export species, uniprotX
export try_download_alphafold, download_alphafolds, getchain
export filter_species!, filter_long!
export project_sequences, columnwise_entropy, align, chargelocations, positive_locations, negative_locations

include("naming_conventions.jl")
include("utils.jl")
include("msa.jl")
include("analyze.jl")

end
