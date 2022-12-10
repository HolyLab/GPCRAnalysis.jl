function chimerax_script(scriptfilename, struct_filenames, ridxs::AbstractVector{<:AbstractVector{<:Integer}};
                         chain_transparency=80, styles=Dict{Tuple{Int,Int},String}(), extras=String[])
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
        if (length(struct_filenames) > 1 &&
            !any(str -> startswith(str, "align"), struct_filenames) &&
            !any(str -> startswith(str, "align"), extras))
            println(io, "matchmaker #2-", length(struct_filenames), " to #1")
        end
        for ex in extras
            println(io, ex)
        end
    end
end

"""
    chimerax_script(scriptfilename, uprot_list, msa::AnnotatedMultipleSequenceAlignment, colidxs;
                    dir=pwd(), chain_transparency=80, styles=Dict{Int,String}(), extras=String[])

Create a [chimerax]() visualization script with name `scriptfilename`. `uprot_list` is a list of UniProtX names that you
want to visualize.  The AlphaFold PDB files for these proteins should be `dir`. `msa` is a Multiple Sequence alignment
and `colidxs` specifies the column indices in `msa` corresponding to amino acid side chains that you'd like to visualize.

`chain_transparency` sets the transparency on the ribbon diagrams (0 = not transparent).
`styles` can be used to affect the display, e.g., `Dict(k => "@SD sphere")` would cause methionines at column index `k`
to be displayed with the sulfur in sphere mode. `extras` can be used to hand-specify a number of additional commands;
this can be useful if, for example, the `msa` has occasional misalignments.

# Examples

Suppose you have the `msa` for rhodopsin (mouse: P15409), then:

```julia
chimerax_script("myscript.cxc", ["P15409"], msa, [i1, i2, i3])
```

where `i1` through `i3` are column-indices in the `msa` that you'd like to view.
"""
function chimerax_script(scriptfilename, uprot_list, msa::AnnotatedMultipleSequenceAlignment, colidxs;
                         dir=pwd(), styles=Dict{Int,String}(), kwargs...)
    ridxs = [Int[] for _ in 1:length(uprot_list)]
    struct_filenames = Vector{String}(undef, length(uprot_list))
    rcstyles = Dict{Tuple{Int,Int},String}()
    for (i, p) in enumerate(uprot_list)
        struct_filenames[i] = joinpath(dir, alphafoldfile(p))
        seqidx = findfirst(str -> startswith(str, p), sequencenames(msa))
        sm = getsequencemapping(msa, seqidx)
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

function markers(modelnum::Integer, positions, radius::Real, color)
    return ["marker #$modelnum position $(pos[1]),$(pos[2]),$(pos[3]) radius $radius color $color" for pos in positions]
end
