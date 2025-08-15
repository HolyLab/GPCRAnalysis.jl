# Biophysical properties of amino acids

"""
    AAProperties

A 3-vector (similar to `SVector{3,Float64}`) that can also be accessed by fields
`charge`, `hydropathy`, and `volume`:
- For standard amino acids, `charge` is -1, 0, or 1, with the exception of
  histidine, which is assigned a charge of 0.1 (assuming pH 7.4).
- `hydropathy` is the Kyte-Doolittle hydropathy index. See Kyte J, Doolittle RF
  (May 1982). "A simple method for displaying the hydropathic character of a
  protein". Journal of Molecular Biology. 157 (1): 105–132.
- `volume` is the van der Waals volume in cubic Angstroms (Å³).
"""
struct AAProperties <: AbstractVector{Float64}
    charge::Float64
    hydropathy::Float64   # from https://www.sciencedirect.com/science/article/pii/0022283682905150?via%3Dihub & https://en.wikipedia.org/wiki/Amino_acid
    volume::Float64       # from http://proteinsandproteomics.org/content/free/tables_1/table08.pdf (van der Waals volume in Å³)
end
Base.size(::AAProperties) = (3,)
Base.getindex(v::AAProperties, i::Int) = i == 1 ? v.charge :
                                         i == 2 ? v.hydropathy :
                                         i == 3 ? v.volume : Base.throw_boundserror(v, i)

StaticArrays.SVector(v::AAProperties) = SVector{3}(v.charge, v.hydropathy, v.volume)

Base.:(+)(v::AAProperties, w::AAProperties) = AAProperties(v.charge + w.charge, v.hydropathy + w.hydropathy, v.volume + w.volume)
Base.:(-)(v::AAProperties, w::AAProperties) = AAProperties(v.charge - w.charge, v.hydropathy - w.hydropathy, v.volume - w.volume)
Base.:(/)(v::AAProperties, r::Real) = AAProperties(v.charge / r, v.hydropathy / r, v.volume / r)
Base.:(*)(u::AAProperties, v::Adjoint{Float64, GPCRAnalysis.AAProperties}) = SVector(u) * SVector(v')'
Base.:(\)(L::LowerTriangular{T, StaticArraysCore.SMatrix{3, 3, T, 9}}, v::AAProperties)  where T<:Real = L \ SVector(v)

const aa_properties = Dict(
    'A' => AAProperties(0, 1.8, 67),
    'R' => AAProperties(1, -4.5, 148),
    'N' => AAProperties(0, -3.5, 96),
    'D' => AAProperties(-1, -3.5, 91),
    'C' => AAProperties(0, 2.5, 86),
    'Q' => AAProperties(0, -3.5, 114),
    'E' => AAProperties(-1, -3.5, 109),
    'G' => AAProperties(0, -0.4, 48),
    'H' => AAProperties(0.1, -3.2, 118),   # at pH 7.4, His has a 10% chance of being charged
    'I' => AAProperties(0, 4.5, 124),
    'L' => AAProperties(0, 3.8, 124),
    'K' => AAProperties(1, -3.9, 135),
    'M' => AAProperties(0, 1.9, 124),
    'F' => AAProperties(0, 2.8, 135),
    'P' => AAProperties(0, -1.6, 90),
    'S' => AAProperties(0, -0.8, 73),
    'T' => AAProperties(0, -0.7, 93),
    'W' => AAProperties(0, -0.9, 163),
    'Y' => AAProperties(0, -1.3, 141),
    'V' => AAProperties(0, 4.2, 105),
)

# z-score the amino acid features using the full covariance matrix
const featμ = sum(p for (r, p) in aa_properties)/length(aa_properties)
const featC = sum((Δp = p - featμ; Δp*Δp') for (r, p) in aa_properties)/(length(aa_properties)-1)
const featChol = cholesky(featC)

const aa_properties_zscored = Dict(r => (featChol.L \ (p - featμ)) for (r, p) in aa_properties)

"""
    P = aa_properties_matrix(msa)

Return a matrix of z-scored biophysical properties (see [`AAProperties`](@ref),
but note that the scaling and interpretation is altered by z-scoring) for each
residue in the MSA. `P` is a matrix-of-vectors.

`P[i, j]` is the property-vector of the `i`th residue (in MSA indexing) of the
`j`th sequence. Note this is transposed relative to the standard MSA matrix.
Transposition facilitates "flattening" the property-vectors along the first axis
using `reinterpret`.
"""
function aa_properties_matrix(msa)
    props = copy(aa_properties_zscored)
    props['-'] = zero(valtype(props))
    props['X'] = zero(valtype(props))
    return [props[Char(residue)] for residue in permutedims(residuematrix(msa))]
end

function aa_properties_matrix(seqs::AbstractVector{FASTX.FASTA.Record})
    props = copy(aa_properties_zscored)
    props['-'] = zero(valtype(props))
    props['X'] = zero(valtype(props))
    return reduce(hcat, [[props[r] for r in sequence(rec)] for rec in seqs])
end