using TEDOPA
using Documenter

DocMeta.setdocmeta!(TEDOPA, :DocTestSetup, :(using TEDOPA); recursive=true)

makedocs(;
    modules=[TEDOPA],
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    sitename="TEDOPA.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
