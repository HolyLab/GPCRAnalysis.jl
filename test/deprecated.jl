@testset "Deprecated" begin
    @info "Testing deprecated functions (warnings are expected)"
    @test only(GPCRAnalysis.markers(2, [Coordinates(5, 4, 3)], 0.5, "blue")) == "marker #2 position 5.0,4.0,3.0 radius 0.5 color blue"
    dot = GPCRAnalysis.marker(2, Coordinates(5, 4, 3), 0.5, "blue")
    tmpfile = tempname() * ".cxc"
    chimerax_script(tmpfile, ["align_ABCD.pdb", "align_EFGH.pdb"], [[5, 10, 15], [7, 11]]; extras=[dot])
    script = read(tmpfile, String)
    @test !occursin("matchmaker #2-2 to #1", script)
    @test occursin("marker #2 position 5.0,4.0,3.0 radius 0.5 color blue", script)
    rm(tmpfile)
end
