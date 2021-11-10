module GPCRAnalysis

using Downloads
using Statistics

using MIToS
using MIToS.Pfam
using MIToS.MSA
using MIToS.PDB

export species, uniprotX
export try_download_alphafold, getchain

include("naming_conventions.jl")
include("utils.jl")

end
