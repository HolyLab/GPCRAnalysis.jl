module GPCRAnalysis

using Downloads
using Statistics
using LinearAlgebra

using MIToS
using MIToS.Pfam
using MIToS.MSA
using MIToS.PDB

using MultivariateStats
using Distances
using Hungarian
using ProgressMeter
using StaticArrays
using Rotations
using TravelingSalesmanHeuristics

using Interpolations
using MutableConvexHulls
using CoordinateTransformations
using GaussianMixtureAlignment

# For querying the Uniprot and Alphafold REST APIs
using HTTP
using JSON
using GZip

export @res_str

export SequenceMapping, AccessionCode, MSACode
export species, uniprotX, query_uniprot_accession
export try_download_alphafold, query_alphafold_latest, download_alphafolds, alphafoldfile, alphafoldfiles, getchain, findall_subseq
export align_to_axes, align_to_membrane
export filter_species!, filter_long!, sortperm_msa, chimerax_script
export project_sequences, columnwise_entropy, align, residue_centroid, residue_centroid_matrix, alphacarbon_coordinates, alphacarbon_coordinates_matrix, mapclosest, chargelocations, positive_locations, negative_locations
export StructAlign, residueindex, ismapped
export BWScheme, lookupbw
export aa_properties, aa_properties_zscored
export sidechaincentroid, scvector, inward_tm_residues, inward_ecl_residues
export features_from_structure

include("naming_conventions.jl")
include("utils.jl")
include("msa.jl")
include("alphafold.jl")
include("analyze.jl")
include("tmalign.jl")
include("align.jl")
include("ballesteros_weinstein.jl")
include("properties.jl")
include("chimerax.jl")
include("pocket.jl")
include("features.jl")

end
