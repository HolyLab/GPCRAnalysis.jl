using GPCRAnalysis
using MIToS
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
    if !isdefined(@__MODULE__, :skip_download) || !skip_download
        @testset "AlphaFold" begin
            @test try_download_alphafold("garbage") === nothing
            fn = tempname()
            @test try_download_alphafold("K7N608", fn) == fn
            @test isfile(fn)
            @test getchain(fn) isa AbstractVector{MIToS.PDB.PDBResidue}
            rm(fn)
        end
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
    end
    @testset "analyze" begin
    end
end
