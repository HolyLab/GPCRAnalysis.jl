struct AAProperties <: AbstractVector{Float64}
    charge::Float64
    hydropathy::Float64   # from https://www.sciencedirect.com/science/article/pii/0022283682905150?via%3Dihub & https://en.wikipedia.org/wiki/Amino_acid
    volume::Float64       # from http://proteinsandproteomics.org/content/free/tables_1/table08.pdf (van der Waals volume in Å³)
end
Base.size(v::AAProperties) = (3,)
Base.getindex(v::AAProperties, i::Int) = i == 1 ? v.charge :
                                         i == 2 ? v.hydropathy :
                                         i == 3 ? v.volume : Base.throw_boundserror(v, i)
const aa_properties = Dict(
    Residue('A') => AAProperties(0, 1.8, 67),
    Residue('R') => AAProperties(1, -4.5, 148),
    Residue('N') => AAProperties(0, -3.5, 96),
    Residue('D') => AAProperties(-1, -3.5, 91),
    Residue('C') => AAProperties(0, 2.5, 86),
    Residue('Q') => AAProperties(0, -3.5, 114),
    Residue('E') => AAProperties(-1, -3.5, 109),
    Residue('G') => AAProperties(0, -0.4, 48),
    Residue('H') => AAProperties(0.1, -3.2, 118),   # at pH 7.4, His has a 10% chance of being charged
    Residue('I') => AAProperties(0, 4.5, 124),
    Residue('L') => AAProperties(0, 3.8, 124),
    Residue('K') => AAProperties(1, -3.9, 135),
    Residue('M') => AAProperties(0, 1.9, 124),
    Residue('F') => AAProperties(0, 2.8, 135),
    Residue('P') => AAProperties(0, -1.6, 90),
    Residue('S') => AAProperties(0, -0.8, 73),
    Residue('T') => AAProperties(0, -0.7, 93),
    Residue('W') => AAProperties(0, -0.9, 163),
    Residue('Y') => AAProperties(0, -1.3, 141),
    Residue('V') => AAProperties(0, 4.2, 105),
)

# z-score the amino acid features using the full covariance matrix
featμ = sum(p for (r, p) in aa_properties)/length(aa_properties)
featC = sum((Δp = p - featμ; Δp*Δp') for (r, p) in aa_properties)/(length(aa_properties)-1)
featChol = cholesky(featC)

const aa_properties_zscored = Dict(r => SVector{3}(featChol.L \ (p - featμ)) for (r, p) in aa_properties)
