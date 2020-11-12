using SDFResultViewer
using Documenter

makedocs(;
    modules=[SDFResultViewer],
    authors="Sebastian Micluța-Câmpeanu <m.c.sebastian95@gmail.com> and contributors",
    repo="https://github.com/SebastianM-C/SDFResultViewer.jl/blob/{commit}{path}#L{line}",
    sitename="SDFResultViewer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SebastianM-C.github.io/SDFResultViewer.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SebastianM-C/SDFResultViewer.jl",
)
