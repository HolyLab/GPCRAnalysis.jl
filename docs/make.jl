using GPCRAnalysis
using Documenter

DocMeta.setdocmeta!(GPCRAnalysis, :DocTestSetup, :(using GPCRAnalysis); recursive=true)

makedocs(;
    modules=[GPCRAnalysis],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    repo="https://github.com/HolyLab/GPCRAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="GPCRAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HolyLab.github.io/GPCRAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HolyLab/GPCRAnalysis.jl",
    devbranch="main",
)
