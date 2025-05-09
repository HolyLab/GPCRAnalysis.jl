using GPCRAnalysis
using MIToS: MSA, Pfam
using MIToS.MSA: @res_str, three2residue, coverage, GappedAlphabet, filtersequences!, percentsimilarity,
                 nsequences, sequencenames, getsequencemapping, getcolumnmapping
using BioStructures
using FASTX
using GaussianMixtureAlignment
using InvertedIndices
using Statistics
using LinearAlgebra
using JuMP, HiGHS
using ColorTypes
using Test

# skip the network-hitting components by setting `skip_download = true` in the global namespace

@testset "GPCRAnalysis.jl" begin
    @test isempty(detect_ambiguities(GPCRAnalysis))
    @testset "utils" begin
        @test @inferred(GPCRAnalysis.residue_range("Q8VGW6_MOUSE/31-308")) == 31:308
        @test @inferred(GPCRAnalysis.strip_residue_range("Q8VGW6_MOUSE/31-308")) == "Q8VGW6_MOUSE"
        @test @inferred(species("Q8VGW6_MOUSE/31-308")) == "MOUSE"
        @test @inferred(uniprotX("Q8VGW6_MOUSE/31-308")) == "Q8VGW6"
        # ClusalOmega MSA codes
        @test uniprotX("Q5QD11|reviewed|Trace") == "Q5QD11"
        # Examples from https://www.uniprot.org/help/accession_numbers
        @test @inferred(uniprotX("A2BC19")) == "A2BC19"
        @test @inferred(uniprotX("P12345")) == "P12345"
        @test @inferred(uniprotX("A0A023GPI8")) == "A0A023GPI8"
        # X identifiers have strict naming conventions (among other things, only lengths 1-5, 6, or 11 are allowed)
        @test_throws ErrorException uniprotX("Q8VGW67_MOUSE/31-308")
        @test_throws ErrorException uniprotX("q8vgw6_MOUSE/31-308")
        @test_throws ErrorException uniprotX("QV8GW6_MOUSE/31-308")
        # FASTA
        @test @inferred(uniprotX("tr|M3WEA8|M3WEA8_FELCA Olfactory receptor family 51 subfamily S member 1 OS=Felis catus OX=9685 GN=OR51S1 PE=4 SV=2")) == "M3WEA8"

        # Versioned
        @test uniprotX("Q8VGW6.37") == "Q8VGW6"
    end

    @testset "NamingConvention" begin
        for NC in (MSACode, AccessionCode)
            @test NC("Q8VGW6") == NC("Q8VGW6")
            @test NC("Q8VGW6") != NC("Q8VGW7")
            @test NC["Q8VGW6"] isa Vector{NC}
            # Disallow implicit conversion, since this loses information
            @test_throws MethodError convert(String, NC("Q8VGW6"))
            # Allow explicit type-casting
            @test String(NC("Q8VGW6")) == "Q8VGW6"
        end
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
        msa = MSA.read_file(pf09645_sto, Pfam.Stockholm)
        @test MSA.nsequences(filter_species!(deepcopy(msa), "ATV")) == 1
        @test MSA.nsequences(filter_long!(deepcopy(msa), 70)) == 3

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
        @test size(project_sequences(msa; fracvar=0.5)) == (1, 4)

        @test AccessionCode(msa, MSACode("Y070_ATV/2-70")) == AccessionCode("Q3V4T1")
        @test MSACode(msa, AccessionCode("Q3V4T1")) == MSACode("Y070_ATV/2-70")
        @test msa[MSACode("Y070_ATV/2-70")][8] == msa[AccessionCode("Q3V4T1")][8] == MSA.Residue('V')
    end
    @testset "Properties" begin
        pf09645_sto = "PF09645_full.stockholm"
        msa = MSA.read_file(pf09645_sto, Pfam.Stockholm)
        X = aa_properties_matrix(msa)
        ΔX = X .- mean(X, dims=2)
        i = findfirst(==(14), getcolumnmapping(msa))
        @test all(iszero, ΔX[i, :])
        @test !all(iszero, ΔX[i-1, :])
        seqs = FASTAReader(open("test.fasta")) do io
            collect(io)
        end
        X = aa_properties_matrix(seqs)
        ΔX = X .- mean(X, dims=2)
        @test all(iszero, ΔX[1, :])
        @test !all(iszero, ΔX[3, :])
    end
    @testset "ChimeraX" begin
        tmpfile = tempname() * ".cxc"
        chimerax_script(tmpfile, ["ABCD.pdb", "EFGH.pdb"], [[5, 10, 15], [7, 11]])
        script = read(tmpfile, String)
        @test occursin("open ABCD", script)
        @test occursin("open EFGH", script)
        @test occursin("show #1 :5", script)
        @test occursin("show #2 :11", script)
        @test occursin("transparency #1-2 80 target c", script)
        @test occursin("matchmaker #2-2 to #1", script)
        dot = GPCRAnalysis.marker(2, [5, 4, 3], 0.5, "blue")
        chimerax_script(tmpfile, ["ABCD.pdb", "EFGH.pdb"], [[5, 10, 15], [7, 11]]; align=false, extras=[dot])
        script = read(tmpfile, String)
        @test !occursin("matchmaker #2-2 to #1", script)
        @test occursin("marker #2 position 5,4,3 radius 0.5 color blue", script)
        rm(tmpfile)
    end
    @testset "Needleman-Wunsch" begin
        function testscore(S, P, D, gapcosts)
            # This tests the premise of dynamic programming, that S[i, j] is the optimal score of align(seq1[1:i], seq2[1:j])
            for i in axes(D, 1), j in axes(D, 2)
                S′, P′ = GPCRAnalysis.score_nw(D[begin:i, begin:j], gapcosts)
                @test S[i, j] == S′[end, end]
                @test P[i, j] == P′[end, end]
            end
        end

        # One-sided matching (must match every element of the first sequence)
        D = [0 1 2 3 4; 1 0 1 2 3; 2 1 0 1 2; 3 2 1 0 1]
        Dboost = D .+ 7 # offset to give something more stringent to test against
        gapcosts = NWGapCosts{Int}(open2=typemax(Int))
        S, P = GPCRAnalysis.score_nw(Dboost, gapcosts)
        # check that disallowed paths are maximally penalized
        for i = 1:4, j = 0:i-1
            @test S[i, j] == typemax(Int)
        end
        @test match(r"4×5 Matrix{GPCRAnalysis.NWParent}:\n ↖  ←  ←  ←  ←\n .  ↖  ←  ←  ←\n .  .  ↖  ←  ←\n .  .  .  ↖  ←"m,
                    sprint(show, MIME("text/plain"), P[1:4, 1:5])) !== nothing
        for ϕ in (GPCRAnalysis.traceback_nw(P), align_nw(D, gapcosts))
            @test ϕ == [(1, 1), (2, 2), (3, 3), (4, 4)]
        end
        ϕ = GPCRAnalysis.traceback_nw(P)
        @test S[end, end] == 28 == sum(Dboost[ij...] for ij in ϕ) + gapcosts(ϕ, axes(D)...)
        testscore(S, P, Dboost, gapcosts)

        D = [0 1 2 3 4; 1 0 1 2 3; 2 1 0 1 2; 4 3 2 1 0]
        @test align_nw(D, gapcosts) == [(1, 1), (2, 2), (3, 3), (4, 5)]
        @test align_nw(D', gapcosts) === nothing   # no way to match all of seq1 if seq1 is longer than seq2
        S, P = GPCRAnalysis.score_nw(D, gapcosts)
        testscore(S, P, D, gapcosts)
        D = [0 1 2 3 4; 1 0 1 2 3; 3 2 1 0 1; 4 3 2 1 0]
        @test align_nw(D, gapcosts) == [(1, 1), (2, 2), (3, 4), (4, 5)]
        D = [1 0 1 2 3; 3 2 1 0 1; 4 3 2 1 0; 5 4 3 2 1]
        @test align_nw(D, gapcosts) == [(1, 2), (2, 3), (3, 4), (4, 5)]

        opsd = getchain("AF-P15409-F1-model_v4.pdb")
        opsd_tms = [37:61, 74:96, 111:133, 153:173, 203:224, 253:274, 287:308]
        @test align_ranges(opsd, opsd, opsd_tms) == opsd_tms

        # Symmetric matching
        D = [1 2 7 4 5;
             2 1 6 3 4;
             3 2 5 2 3;]
        gapcosts = NWGapCosts{Int}(extend1=1, extend2=1, open1=3, open2=3)
        S, P = GPCRAnalysis.score_nw(D, gapcosts)
        ϕ = GPCRAnalysis.traceback_nw(P)
        @test ϕ == [(1, 1), (2, 2), (3, 5)]
        @test S[end, end] == sum(D[ij...] for ij in ϕ) + gapcosts(ϕ, axes(D)...)
        testscore(S, P, D, gapcosts)
    end
    @testset "Pocket residues and features" begin
        opsd = getchain("AF-P15409-F1-model_v4.pdb")
        opsdr = [three2residue(String(resname(r))) for r in opsd]
        # These are not quite accurate, from https://www.uniprot.org/uniprotkb/P15409/entry#subcellular_location
        # the actual answer is
        #   opsd_tms = [37:61, 74:96, 111:133, 153:173, 203:224, 253:274, 287:308]
        # but these exercise useful parts of the code
        opsd_tms = [only(findall_subseq(res"PWQF", opsdr)):only(findall_subseq(res"VTVQ", opsdr))+3,
                    only(findall_subseq(res"NYIL", opsdr)):only(findall_subseq(res"YTSL", opsdr))+3,
                    only(findall_subseq(res"PTGC", opsdr)):only(findall_subseq(res"YVVV", opsdr))+3,
                    only(findall_subseq(res"ENHA", opsdr)):only(findall_subseq(res"PPLV", opsdr))+3,
                    only(findall_subseq(res"NESF", opsdr)):only(findall_subseq(res"LVFT", opsdr))+3,
                    only(findall_subseq(res"AEKE", opsdr)):only(findall_subseq(res"YIFT", opsdr))+3,
                    only(findall_subseq(res"PIFM", opsdr)):only(findall_subseq(res"YIML", opsdr))+3,
        ]
        applytransform!(opsd, align_to_membrane(opsd, opsd_tms))
        # Check the alignment
        leaflet_e, leaflet_i = GPCRAnalysis.collect_leaflet_residues(opsd, opsd_tms, true)
        ce = reduce(hcat, [coordarray(r) for r in leaflet_e])
        ci = reduce(hcat, [coordarray(r) for r in leaflet_i])
        ze = median(ce[3,:])
        zi = median(ci[3,:])
        @test ze > zi    # extracellular is positive
        @test abs(ze + zi) < 1e-6 * (ze - zi) # center of membrane is at z=0

        tm_res = inward_tm_residues(opsd, opsd_tms[[2,3,5,6,7]])
        @test length(tm_res) == 5
        for (i,tm) in enumerate(tm_res)
            @test length(tm) == length(opsd_tms[[2,3,5,6,7]][i])
            @test isa(tm, Vector{Bool})
        end

        opsd_ecls = [1:opsd_tms[1][1]-1,
                     opsd_tms[2][end]+1:opsd_tms[3][1]-1,
                     opsd_tms[4][end]+1:opsd_tms[5][1]-1,
                     opsd_tms[6][end]+1:opsd_tms[7][1]-1,
        ]
        ecl_res = inward_ecl_residues(opsd, opsd_ecls)
        @test length(ecl_res) == length(opsd_ecls)
        for (i,ecl) in enumerate(ecl_res)
            @test length(ecl) == length(opsd_ecls[i])
            @test isa(ecl, Vector{Bool})
        end

        tm_idxs = vcat([opsd_tms[[2,3,5,6,7]][i][tm_res[i]] for i=1:5]...)
        tm_mgmm = features_from_structure(opsd, tm_idxs)
        tm_mgmm_combined = features_from_structure(opsd, tm_idxs; combined=true)
        tm_mgmm_combined_2 = features_from_structure(opsd, tm_idxs; combined=true, σfun=GPCRAnalysis.equalvolumeradius, ϕfun=(a,r,f)->1.0)
        @test sum(length, values(tm_mgmm.gmms)) > sum(length, values(tm_mgmm_combined.gmms))
        for key in keys(tm_mgmm.gmms)
            for (g1, g2) in zip(tm_mgmm_combined.gmms[key].gaussians, tm_mgmm_combined_2.gmms[key].gaussians)
                @test g1.μ == g2.μ
                @test g1.σ >= g2.σ
                @test g1.ϕ >= g2.ϕ
            end
        end

        ecl_idxs = vcat([opsd_ecls[i][ecl_res[i]] for i=1:4]...)
        ecl_mgmm = features_from_structure(opsd, ecl_idxs)
        ecl_mgmm_combined = features_from_structure(opsd, ecl_idxs; combined=true)
        ecl_mgmm_combined_2 = features_from_structure(opsd, ecl_idxs; combined=true, σfun=GPCRAnalysis.equalvolumeradius, ϕfun=(a,r,f)->1.0)
        @test sum(length, values(ecl_mgmm.gmms)) > sum(length, values(ecl_mgmm_combined.gmms))
        for key in keys(ecl_mgmm.gmms)
            for (g1, g2) in zip(ecl_mgmm_combined.gmms[key].gaussians, ecl_mgmm_combined_2.gmms[key].gaussians)
                @test g1.μ == g2.μ
                @test g1.σ >= g2.σ
                @test g1.ϕ >= g2.ϕ
            end
        end
    end

    @testset "Forces" begin
        interactions = [(:Steric, :Steric) => 1, (:Hydrophobe, :Hydrophobe) => -1]
        mgmmx = IsotropicMultiGMM(Dict(:Steric     => IsotropicGMM([IsotropicGaussian([0.0, 0.0], 1.0, 1.0)]),
                                       :Hydrophobe => IsotropicGMM([IsotropicGaussian([0.0, 0.0], 2.0, 1.0)])))
        mgmmy = IsotropicMultiGMM(Dict(:Steric     => IsotropicGMM([IsotropicGaussian([1.0, 1.0], 1.0, 1.0)]),
                                       :Hydrophobe => IsotropicGMM([IsotropicGaussian([1.0, 1.0], 2.0, 1.0)])))
        forces = forcecomponents([mgmmx, mgmmy], interactions)
        @test length(forces) == 2
        @test iszero(forces[1] + forces[2])
        @test all(f -> norm(sum(f; dims=2)) > 0.1, forces)   # non-zero force
        w = optimize_weights(forces)
        @test all(>=(0), w)
        @test sum(w) ≈ 1 atol=1e-6
        badinteractions = [(:Steric, :Steric) => 1, (:Hydrophobe, :Hydrophobe) => -1, (:Ionic, :Ionic) => -1]
        @test_throws "is not sign-compatible" forcecomponents([mgmmx, mgmmy], badinteractions)
        weighted_interactions = [key => w0*wi for ((key, w0), wi) in zip(interactions, w)]
        wforces = forcecomponents([mgmmx, mgmmy], weighted_interactions)
        @test all(f .* w' ≈ wf for (f, wf) in zip(forces, wforces))
        @test all(f -> norm(sum(f; dims=2)) < 1e-4, wforces)   # net-zero force
        wnew = optimize_weights(wforces)
        @test wnew[1] ≈ 0.5 atol=1e-4
        @test wnew[2] ≈ 0.5 atol=1e-4

        #
        #      Y
        #          W   X
        #      Z
        #
        #where WX is hydrophobic, WY is H-bond, and WZ is ionic (all attractive)
        interactions = [(:Hydrophobe, :Hydrophobe) => -1, (:Donor, :Acceptor) => -1, (:Ionic, :Ionic) => 1]
        mgmmw = IsotropicMultiGMM(Dict(:Hydrophobe   => IsotropicGMM([IsotropicGaussian([ 0.0, 0.0], 2.0, 1.0)]),
                                       :Donor        => IsotropicGMM([IsotropicGaussian([ 0.0, 0.0], 1.5, 1.0)]),
                                       :PosIonizable => IsotropicGMM([IsotropicGaussian([ 0.0, 0.0], 4.0, 1.0)])))
        mgmmx = IsotropicMultiGMM(Dict(:Hydrophobe   => IsotropicGMM([IsotropicGaussian([ 1.0, 0.0], 2.0, 1.0)])))
        mgmmy = IsotropicMultiGMM(Dict(:Acceptor     => IsotropicGMM([IsotropicGaussian([-1.0, 1.0], 1.5, 1.0)])))
        mgmmz = IsotropicMultiGMM(Dict(:NegIonizable => IsotropicGMM([IsotropicGaussian([-1.0,-1.0], 4.0, 1.0)])))
        forces = forcecomponents([mgmmw, mgmmx, mgmmy, mgmmz], interactions, [1])
        @test length(forces) == 1
        @test all(f -> norm(sum(f; dims=2)) > 0.1, forces)   # non-zero force
        f = only(forces)
        F = f' * f
        # Test that each force is in the expected direction
        @test f[:, 1] ≈ ([ 1, 0] \ f[:, 1]) * [ 1, 0]
        @test f[:, 2] ≈ ([-1, 1] \ f[:, 2]) * [-1, 1]
        @test f[:, 3] ≈ ([-1,-1] \ f[:, 3]) * [-1,-1]
        w = optimize_weights(forces)
        @test all(>=(0), w)
        @test sum(w) ≈ 1
        weighted_interactions = [key => w0*wi for ((key, w0), wi) in zip(interactions, w)]
        wforces = forcecomponents([mgmmw, mgmmx, mgmmy, mgmmz], weighted_interactions, [1])
        @test all(f .* w' ≈ wf for (f, wf) in zip(forces, wforces))
        @test all(f -> norm(sum(f; dims=2)) < 1e-6, wforces)   # net-zero force
        wF = sum(f' * f for f in wforces)
        @test all(abs.(wF*ones(3)) .< 1e-7)   # gradient is zero at the minimum even under re-tuning

        interactions = [(:Steric, :Steric) => 1, (:Hydrophobe, :Hydrophobe) => -1, (:Donor, :Acceptor) => -1]  # drop the ionic
        opsd = getchain("AF-P15409-F1-model_v4.pdb")
        opsd_tms = [37:61, 74:96, 111:133, 153:173, 203:224, 253:274, 287:308]
        tm_res = inward_tm_residues(opsd, opsd_tms)
        tm_idxs = reduce(vcat, [tm[idx] for (tm, idx) in zip(opsd_tms, tm_res)])
        for combined in (false, true)
            forces = forcecomponents(opsd, interactions, tm_idxs; combined)
            w = optimize_weights(forces)
            # Exact force balance is presumably impossible to achieve in this case.
            # For a weaker test, check that weighted interactions are scaled as expected.
            weighted_interactions = [key => w0*wi for ((key, w0), wi) in zip(interactions, w)]
            wforces = forcecomponents(opsd, weighted_interactions, tm_idxs; combined)
            @test all(f .* w' ≈ wf for (f, wf) in zip(forces, wforces))
        end
    end

    if !isdefined(@__MODULE__, :skip_download) || !skip_download
        @testset "Uniprot" begin
            @test query_uniprot_accession("T2R38_MOUSE") == "Q7TQA6"

            jobID = GPCRAnalysis.map_uniprot_submit(["T2R38_MOUSE"])
            tstart = time()
            sleep(1)
            status = false
            while status != true && time() < tstart + 20
                sleep(1)
                status = GPCRAnalysis.map_uniprot_status(jobID)
            end
            results = GPCRAnalysis.map_uniprot_retrieve(jobID)[:results]
            resultsdict = Dict(obj["from"] => obj["to"] for obj in results)
            @test resultsdict["T2R38_MOUSE"] == "Q7TQA6"

            jobID = GPCRAnalysis.map_uniprot_submit(["ENSMUSG00000067064", "ENSMUSG00000057464"], "Ensembl")
            tstart = time()
            sleep(1)
            status = false
            while status != true && time() < tstart + 20
                sleep(1)
                status = GPCRAnalysis.map_uniprot_status(jobID)
            end
            results = GPCRAnalysis.map_uniprot_retrieve(jobID)[:results]
            resultsdict = Dict(obj["from"] => obj["to"] for obj in results)
            @test resultsdict["ENSMUSG00000067064"] == "Q8VGU4"


        end

        @testset "EBI and NCBI" begin
            result = only(query_ebi_proteins("Q7TQA6"))
            obj = nothing
            for x in result["dbReferences"]
                if x["type"] == "Proteomes"
                    obj = x
                    break
                end
            end
            m = match(r"Chromosome\s(\d+)", obj["properties"]["component"])
            @test parse(Int, m.captures[1]) == 6
            for x in result["dbReferences"]
                if x["type"] == "RefSeq"
                    obj = x["id"]
                    break
                end
            end
            q = query_ncbi(obj)
            @test q["total_count"] == 1
        end

        @testset "AlphaFold" begin
            @test try_download_alphafold("garbage") === nothing
            uname = "K7N608"
            url = query_alphafold_latest(uname)
            v = parse(Int, match(GPCRAnalysis.rex_alphafold_pdbs, url).captures[2])
            mktempdir() do dir
                fn = split(url, '/')[end]
                absfn = joinpath(dir, fn)
                @test try_download_alphafold(uname, absfn; version=v) == absfn
                @test isfile(absfn)
                @test getchain(absfn) isa GPCRAnalysis.ChainLike
                @test only(alphafoldfiles(dir)) == fn
            end
        end

        @testset "MSA + AlphaFold" begin
            mktempdir() do path
                pfamfile = "PF03402.alignment.full.gz"
                pfampath = joinpath(path, pfamfile)
                if !isfile(pfampath)
                    Pfam.downloadpfam("PF03402"; filename=pfampath)
                end
                msa = MSA.read_file(pfampath, MSA.Stockholm, generatemapping=true, useidcoordinates=true)
                filter_species!(msa, "MOUSE")
                # Make small enough for a decent download test
                filtersequences!(msa, coverage(msa) .>= 0.9)
                filtersequences!(msa, startswith.(names(msa.matrix, 1), Ref("K7N7")))

                sim0 = percentsimilarity(msa); sim0 = sum(sim0[i-1, i] for i = 2:nsequences(msa))
                tour = sortperm_msa(msa)
                @test isperm(tour)
                msa = msa[tour,:]
                sim = percentsimilarity(msa); sim = sum(sim[i-1, i] for i = 2:nsequences(msa))
                @test sim >= sim0

                msaentropies = columnwise_entropy(msa)
                conserved_cols = findall(msaentropies .< 0.5)

                download_alphafolds(msa; dirname=path)
                nfiles = length(readdir(path))
                download_alphafolds(["K7N701", "G3UY47"]; dirname=path)
                @test length(readdir(path)) == nfiles + 1
                msacode2structfile = alphafoldfiles(msa, path)
                afnbyidx(i) = getchain(joinpath(path, msacode2structfile[MSACode(sequencenames(msa)[i])]))

                c1, c2 = afnbyidx(1), afnbyidx(5)
                @test @inferred(pLDDTcolor(first(c1))) isa AbstractRGB
                conserved_residues = c1[SequenceMapping(getsequencemapping(msa, 1))[conserved_cols]]
                badidx = findall(==(nothing), conserved_residues)
                conserved_residues = convert(Vector{Residue}, conserved_residues[Not(badidx)])
                conserved_cols = conserved_cols[Not(badidx)]
                sm = SequenceMapping(getsequencemapping(msa, 2))
                c2a = applytransform!(copy(c2), align(conserved_residues, c2, sm[conserved_cols]))
                @test rmsd(conserved_residues, skipnothing(c2a[sm[conserved_cols]]); superimpose=false) <
                      rmsd(conserved_residues, skipnothing(c2[sm[conserved_cols]]); superimpose=false)
                mc = map_closest(c1, c2a)
                @test align_closest(c1, c2a) isa Transformation

                cloc = chargelocations(c1)
                lpos = positive_locations(cloc)
                @test !isempty(lpos) && all(item -> isa(item, AbstractVector), lpos)
                lneg = negative_locations(cloc)
                @test !isempty(lneg) && all(item -> isa(item, AbstractVector), lneg)
                @test isempty(lpos ∩ lneg)

                # Choose a sufficiently-divergent pair that structural alignment is nontrivial
                idxref = findfirst(str -> startswith(str, "K7N701"), sequencenames(msa))
                idxcmp = findfirst(str -> startswith(str, "K7N778"), sequencenames(msa))
                cref, ccmp = afnbyidx(idxref), afnbyidx(idxcmp)
                sa = StructAlign(cref, ccmp, joinpath(@__DIR__, "tmalign.txt"))
                @test !ismapped(sa, 1, nothing)
                @test  ismapped(sa, 11, nothing)
                @test  ismapped(sa, nothing, 1)
                @test !ismapped(sa, nothing, length(ccmp))
                @test_throws BoundsError ismapped(sa, nothing, length(ccmp)+1)
                @test_throws BoundsError residueindex(sa, 1, nothing)
                @test residueindex(sa, 2, nothing, 1) == 1
                @test residueindex(sa, 12, nothing, -1) == 1
                @test residueindex(sa, nothing, 1) == 11
                @test residueindex(sa, nothing, length(ccmp), -1) == 304
                @test residueindex(sa, 304, nothing) == 296
                conserved_cols = findall(msaentropies .< 0.5)
                smref = SequenceMapping(getsequencemapping(msa, idxref))[conserved_cols]
                keep = (!iszero).(smref)
                conserved_cols = conserved_cols[keep]
                smref = smref[keep]
                conserved_residues = cref[smref]
                smcmp = SequenceMapping(getsequencemapping(msa, idxcmp))
                ccmpa = applytransform!(copy(ccmp), align(conserved_residues, ccmp, smcmp[conserved_cols]))
                mc = map_closest(cref, ccmpa)
                idxclose = first.(filter(item -> item[2] < 5.0, mc))   # TMAlign scores those closer than 5Å as a match
                n = 0; for i in Iterators.drop(eachindex(idxclose), 1)
                     n += idxclose[i] < idxclose[i-1]
                end
                @test n < 0.05 * length(idxclose)  # almost sorted
                @test length(setdiff(idxclose, sa.m2.a2s)) < 0.05 * length(idxclose)  # almost same as TMalign

                chimerafile = tempname() * ".cxc"
                chimerax_script(chimerafile, ["K7N775", "K7N731"], msa, [66, 69]; dir=path)
                script = read(chimerafile, String)
                @test  occursin("open $(joinpath(path, "AF-K7N775"))", script)
                @test  occursin("open $(joinpath(path, "AF-K7N731"))", script)
                @test !occursin("open $(joinpath(path, "AF-K7N701"))", script)
                @test occursin("show #1 :67", script)
                @test occursin("show #1 :70", script)
                @test occursin("show #2 :70", script)
                @test occursin("show #2 :73", script)
                @test occursin("matchmaker #2-2 to #1", script)
            end
        end
    end
end
