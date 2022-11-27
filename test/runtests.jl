using GPCRAnalysis
using MIToS
using MIToS.MSA
using MIToS.PDB
using InvertedIndices
using Test

# skip the network-hitting components by setting `skip_download = true` in the global namespace

@testset "GPCRAnalysis.jl" begin
    @test isempty(detect_ambiguities(GPCRAnalysis))
    @testset "utils" begin
        @test @inferred(GPCRAnalysis.residue_range("Q8VGW6_MOUSE/31-308")) == 31:308
        @test @inferred(GPCRAnalysis.strip_residue_range("Q8VGW6_MOUSE/31-308")) == "Q8VGW6_MOUSE"
        @test @inferred(species("Q8VGW6_MOUSE/31-308")) == "MOUSE"
        @test @inferred(uniprotX("Q8VGW6_MOUSE/31-308")) == "Q8VGW6"
        # Examples from https://www.uniprot.org/help/accession_numbers
        @test @inferred(uniprotX("A2BC19")) == "A2BC19"
        @test @inferred(uniprotX("P12345")) == "P12345"
        @test @inferred(uniprotX("A0A023GPI8")) == "A0A023GPI8"
        # X identifiers have strict naming conventions (among other things, only lengths 1-5, 6, or 11 are allowed)
        @test_throws ErrorException uniprotX("Q8VGW67_MOUSE/31-308")
        @test_throws ErrorException uniprotX("q8vgw6_MOUSE/31-308")
        @test_throws ErrorException uniprotX("QV8GW6_MOUSE/31-308")
    end

    @testset "Ballesteros-Weinstein" begin
        opsd_scheme = BWScheme([55, 83, 135, 161, 215, 267, 303],
                               [34:64, 73:99, 107:139, 150:173, 200:229, 246:277, 285:309])
        @test lookupbw(160, opsd_scheme) == (4, 49)
        @test lookupbw((4, 49), opsd_scheme) == 160
    end

    @testset "MSA" begin
        # The test file is copied from MIToS/test/data, with gratitude
        pf09645_sto = "PF09645_full.stockholm"
        msa = read(pf09645_sto, MIToS.Pfam.Stockholm)
        @test MIToS.MSA.nsequences(filter_species!(deepcopy(msa), "ATV")) == 1
        @test MIToS.MSA.nsequences(filter_long!(deepcopy(msa), 70)) == 3

        idx = SequenceMapping([0, 4, 5, 0])
        seqvals = fill(NaN, 9)
        seqvals[idx] = [0.1, 0.2, 0.3, 0.4]
        @test seqvals[4] == 0.2
        @test seqvals[5] == 0.3
        @test all(isnan, seqvals[1:3])
        @test all(isnan, seqvals[6:end])

        # analyze
        e = columnwise_entropy(msa)
        @test length(e) == size(msa, 2) && e[9] == 0
        e2 = columnwise_entropy(msa, GappedAlphabet())
        @test all(e2 .>= e)
        @test !all(e2 .== e)

        @test size(project_sequences(msa)) == (3, 4)
    end

    if !isdefined(@__MODULE__, :skip_download) || !skip_download
        @testset "AlphaFold" begin
            @test try_download_alphafold("garbage") === nothing
            fn = tempname()
            @test try_download_alphafold("K7N608", fn) == fn
            @test isfile(fn)
            @test getchain(fn) isa AbstractVector{MIToS.PDB.PDBResidue}
            rm(fn)
            mktempdir() do path
                pfamfile = "PF03402.alignment.full.gz"
                pfampath = joinpath(path, pfamfile)
                if !isfile(pfampath)
                    Pfam.downloadpfam("PF03402"; filename=pfampath)
                end
                msa = read(pfampath, MSA.Stockholm, generatemapping=true, useidcoordinates=true)
                filter_species!(msa, "MOUSE")
                # Make small enough for a decent download test
                filtersequences!(msa, coverage(msa) .>= 0.9)
                filtersequences!(msa, startswith.(names(msa.matrix, 1), Ref("K7N7")))

                msaentropies = columnwise_entropy(msa)
                conserved_cols = findall(msaentropies .< 0.5)

                download_alphafolds(msa; dirname=path)
                fns = filter(fn -> endswith(fn, ".pdb"), readdir(path; join=true))
                @test length(fns) == nsequences(msa)
                p = sortperm(sequencenames(msa))
                fns = fns[invperm(p)]  # now the filenames correspond to the rows of msa

                c1, c5 = getchain(fns[1]), getchain(fns[5])
                conserved_residues = c1[SequenceMapping(getsequencemapping(msa, 1))[conserved_cols]]
                badidx = findall(==(nothing), conserved_residues)
                conserved_residues = convert(Vector{PDBResidue}, conserved_residues[Not(badidx)])
                conserved_cols = conserved_cols[Not(badidx)]
                sm = SequenceMapping(getsequencemapping(msa, 5))
                c5a = align(conserved_residues, c5, sm[conserved_cols])
                @test rmsd(conserved_residues, c5a[sm[conserved_cols]]; superimposed=true) < rmsd(conserved_residues, c5[sm[conserved_cols]]; superimposed=true)

                cl = chargelocations(c1)
                lpos = positive_locations(cl)
                @test !isempty(lpos) && all(item -> isa(item, Coordinates), lpos)
                lneg = negative_locations(cl)
                @test !isempty(lneg) && all(item -> isa(item, Coordinates), lneg)
                @test isempty(lpos ∩ lneg)

                c3 = getchain(fns[3])
                sa = StructAlign(c3, c5, joinpath(@__DIR__, "tmalign_3_5.txt"))
                @test !ismapped(sa, 1, nothing)
                @test  ismapped(sa, 11, nothing)
                @test  ismapped(sa, nothing, 1)
                @test !ismapped(sa, nothing, length(c5))
                @test_throws BoundsError ismapped(sa, nothing, length(c5)+1)
                @test_throws BoundsError residueindex(sa, 1, nothing)
                @test residueindex(sa, 2, nothing, 1) == 1
                @test residueindex(sa, 12, nothing, -1) == 1
                @test residueindex(sa, nothing, 1) == 11
                @test residueindex(sa, nothing, length(c5), -1) == 304
                @test residueindex(sa, 304, nothing) == 296
            end
        end
    end
end
