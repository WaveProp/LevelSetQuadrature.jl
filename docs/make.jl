using LevelSetQuadrature
using Documenter

DocMeta.setdocmeta!(LevelSetQuadrature, :DocTestSetup, :(using LevelSetQuadrature); recursive=true)

makedocs(;
    modules=[LevelSetQuadrature],
    authors="Luiz M. Faria <maltezfaria@gmail.com> and contributors",
    repo="https://github.com/WaveProp/LevelSetQuadrature.jl/blob/{commit}{path}#{line}",
    sitename="LevelSetQuadrature.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://WaveProp.github.io/LevelSetQuadrature.jl",
        assets=String[],
    ),
    pages=[
        "Getting started" => "index.md",
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/WaveProp/LevelSetQuadrature.jl",
    devbranch="main"
)
