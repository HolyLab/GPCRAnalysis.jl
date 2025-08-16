module GPCRAnalysis

using Downloads
using Statistics
using LinearAlgebra

using BioStructures
using BioStructures: amino_acid_data
using FASTX

using MultivariateStats
using Distances
using OffsetArrays
using Hungarian
using ProgressMeter
using StaticArrays
using Rotations
using TravelingSalesmanHeuristics

using Interpolations
using MutableConvexHulls
using GaussianMixtureAlignment

# For querying the Uniprot and Alphafold REST APIs
using HTTP
using JSON3
using GZip

# AlphaFold pLDDT colors
using FixedPointNumbers
using ColorTypes

const ChainLike = Union{Chain, AbstractVector{<:AbstractResidue}}   # an iterable of residues
const ResidueLike = Union{Residue, AbstractVector{<:AbstractAtom}}  # an iterable of atoms
const StructureLike = Union{ChainLike, Model, MolecularStructure}

# export @res_str

export SequenceMapping, AccessionCode, MSACode, NWGapCosts
export sequenceindexes, columnindexes, isgap, isunknown, sequencekeys, msasequence, residuematrix, subseqs, subseqs!
export species, uniprotX, query_uniprot_accession, query_ebi_proteins, query_ncbi
export try_download_alphafold, query_alphafold_latest, download_alphafolds, alphafoldfile, alphafoldfiles, getchain,
       findall_subseq, pLDDT, pLDDTcolor
export align_to_axes, align_to_membrane, align_nw, align_ranges, map_closest, align_closest
export filter_species!, filter_long!, sortperm_msa, chimerax_script
export project_sequences, columnwise_entropy, align, residue_centroid, residue_centroid_matrix, alphacarbon_coordinates,
       alphacarbon_coordinates_matrix, chargelocations, positive_locations, negative_locations
export StructAlign, residueindex, ismapped, skipnothing
export BWScheme, lookupbw
export aa_properties, aa_properties_zscored, aa_properties_matrix
export sidechaincentroid, scvector, inward_tm_residues, inward_ecl_residues
export features_from_structure
export forcecomponents, optimize_weights, forcedict

include("consts.jl")
include("naming_conventions.jl")
include("msa.jl")
include("utils.jl")
include("query.jl")
include("alphafold.jl")
include("analyze.jl")
include("tmalign.jl")
include("align.jl")
include("ballesteros_weinstein.jl")
include("properties.jl")
include("chimerax.jl")
include("pocket.jl")
include("features.jl")
include("forces.jl")

end
