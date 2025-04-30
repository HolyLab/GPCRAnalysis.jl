# Biophysical properties of amino acids

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
    MSA.Residue('A') => AAProperties(0, 1.8, 67),
    MSA.Residue('R') => AAProperties(1, -4.5, 148),
    MSA.Residue('N') => AAProperties(0, -3.5, 96),
    MSA.Residue('D') => AAProperties(-1, -3.5, 91),
    MSA.Residue('C') => AAProperties(0, 2.5, 86),
    MSA.Residue('Q') => AAProperties(0, -3.5, 114),
    MSA.Residue('E') => AAProperties(-1, -3.5, 109),
    MSA.Residue('G') => AAProperties(0, -0.4, 48),
    MSA.Residue('H') => AAProperties(0.1, -3.2, 118),   # at pH 7.4, His has a 10% chance of being charged
    MSA.Residue('I') => AAProperties(0, 4.5, 124),
    MSA.Residue('L') => AAProperties(0, 3.8, 124),
    MSA.Residue('K') => AAProperties(1, -3.9, 135),
    MSA.Residue('M') => AAProperties(0, 1.9, 124),
    MSA.Residue('F') => AAProperties(0, 2.8, 135),
    MSA.Residue('P') => AAProperties(0, -1.6, 90),
    MSA.Residue('S') => AAProperties(0, -0.8, 73),
    MSA.Residue('T') => AAProperties(0, -0.7, 93),
    MSA.Residue('W') => AAProperties(0, -0.9, 163),
    MSA.Residue('Y') => AAProperties(0, -1.3, 141),
    MSA.Residue('V') => AAProperties(0, 4.2, 105),
)

# z-score the amino acid features using the full covariance matrix
const featμ = sum(p for (r, p) in aa_properties)/length(aa_properties)
const featC = sum((Δp = p - featμ; Δp*Δp') for (r, p) in aa_properties)/(length(aa_properties)-1)
const featChol = cholesky(featC)

const aa_properties_zscored = Dict(r => (featChol.L \ (p - featμ)) for (r, p) in aa_properties)

function aa_properties_matrix(msa::AbstractMultipleSequenceAlignment)
    props = copy(aa_properties_zscored)
    props[GAP] = zero(valtype(props))
    props[MSA.Residue('X')] = zero(valtype(props))
    return [props[residue] for residue in permutedims(msa)]
end

function aa_properties_matrix(seqs::AbstractVector{FASTX.FASTA.Record})
    props = Dict(Char(r) => v for (r, v) in aa_properties_zscored)
    props['-'] = zero(valtype(props))
    props['X'] = zero(valtype(props))
    return reduce(hcat, [[props[r] for r in sequence(rec)] for rec in seqs])
end
