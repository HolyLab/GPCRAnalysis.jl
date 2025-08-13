function chimerax_script(scriptfilename, struct_filenames, ridxs::AbstractVector{<:AbstractVector{<:Integer}};
                         align::Bool=true, chain_transparency=80, styles=Dict{Tuple{Int,Int},String}(), extras=String[])
    open(scriptfilename, "w") do io
        for (i, afname) in enumerate(struct_filenames)
            println(io, "open ", afname)
            for ridx in ridxs[i]
                ridx == 0 && continue
                println(io, "show #", i, " :", ridx)
                style = get(styles, (i, ridx), nothing)
                if style !== nothing
                    println(io, "style #", i, " :", ridx, " ", style)
                end
            end
        end
        println(io, "transparency #1-", length(struct_filenames), ' ', chain_transparency, " target c")
        if (length(struct_filenames) > 1 && align)
            if (any(str -> startswith(str, "align"), struct_filenames) ||
                any(str -> startswith(str, "align"), extras))
                @warn "skipping alignment based on file names is deprecated, use `align=false` instead" maxlog=1
            else
                println(io, "matchmaker #2-", length(struct_filenames), " to #1")
            end
        end
        if isa(extras, AbstractString)
            println(io, extras)
        else
            for ex in extras
                println(io, ex)
            end
        end
    end
end

"""
    chimerax_script(scriptfilename, uprot_list, msa, colidxs;
                    dir=pwd(), align=true, chain_transparency=80, styles=Dict{Int,String}(), extras=String[])

Create a [chimerax](https://www.cgl.ucsf.edu/chimerax/) visualization script
with name `scriptfilename`. `uprot_list` is a list of UniProtX names that you
want to visualize. `msa` is a Multiple Sequence alignment and `colidxs`
specifies the column indices in `msa` corresponding to amino acid side chains
that you'd like to visualize.

Keyword arguments:
- `dir` is the directory with the protein structure files
- `align` determines whether to align the structures to the first one (uses the
  `matchmaker` tool)
- `chain_transparency` sets the transparency on the ribbon diagrams (0 = not
transparent)
- `styles` can be used to affect the display, e.g., `Dict(k => "@SD sphere")`
would cause methionines at column index `k` to be displayed with the sulfur in
sphere mode.
- `extras` can be used to hand-specify a number of additional commands; this can
be useful if, for example, the `msa` has occasional misalignments.

# Examples

Suppose you have the `msa` for rhodopsin (mouse: P15409), then:

```julia
chimerax_script("myscript.cxc", ["P15409"], msa, [i1, i2, i3])
```

where `i1` through `i3` are column-indices in the `msa` that you'd like to view.
"""
function chimerax_script(scriptfilename, uprot_list, msa, colidxs;
                         dir=pwd(), styles=Dict{Int,String}(), kwargs...)
    ridxs = [Int[] for _ in 1:length(uprot_list)]
    struct_filenames = Vector{String}(undef, length(uprot_list))
    rcstyles = Dict{Tuple{Int,Int},String}()
    afs = alphafoldfiles(msa, dir; join=true)
    uprot2msaidx = Dict{AccessionCode,Int}(AccessionCode(msa, name) => i for (i, name) in enumerate(sequencekeys(msa)))
    for (i, p) in enumerate(uprot_list)
        j = uprot2msaidx[AccessionCode(p)]
        struct_filenames[i] = afs[MSACode(sequencekeys(msa)[j])]
        sm = sequenceindexes(msa, j)
        for (j, c) in enumerate(colidxs)
            ridx = sm[c]
            if iszero(ridx)
                @warn "column $c not set in $p"
                continue
            end
            push!(ridxs[i], ridx)
            style = get(styles, c, nothing)
            if style !== nothing
                rcstyles[(i, j)] = style
            end
        end
    end
    return chimerax_script(scriptfilename, struct_filenames, ridxs; styles=rcstyles, kwargs...)
end

"""
    str = marker(modelnum, pos, radius, color)

Return a string that represents a marker in ChimeraX. `modelnum` is the model number, `pos` is a 3D position,
`radius` is the radius of the marker, and `color` is the color of the marker.
"""
@noinline function marker(modelnum::Integer, pos::AbstractVector{<:Real}, radius::Real, color)
    return "marker #$modelnum position $(pos[1]),$(pos[2]),$(pos[3]) radius $radius color $color"
end
