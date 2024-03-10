using TEDOPA
using Documenter

DocMeta.setdocmeta!(TEDOPA, :DocTestSetup, :(using TEDOPA); recursive=true)

makedocs(;
    modules=[TEDOPA],
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    sitename="TEDOPA.jl",
    repo=Remotes.GitHub("phaerrax", "TEDOPA.jl"),
    format=Documenter.HTML(; edit_link="github"),
    pages=["Home" => "index.md"],
    checkdocs=:none,
    deploydocs(; repo="github.com/phaerrax/TEDOPA.jl.git"),
)
