module GPCRAnalysis

using Downloads
using Statistics

using MIToS
using MIToS.Pfam
using MIToS.MSA
using MIToS.PDB

export species, uniprotX
export try_download_alphafold, getchain
export filter_species!, filter_long!

include("naming_conventions.jl")
include("utils.jl")
include("msa.jl")

end
