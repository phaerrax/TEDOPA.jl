using Documenter, DocumenterCitations
using TEDOPA

push!(LOAD_PATH, "../src/")
# This line is needed if the package is not accessible through Julia's LOAD_PATH.

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(;
    modules=[TEDOPA],
    sitename="TEDOPA.jl",
    repo=Remotes.GitHub("phaerrax", "TEDOPA.jl"),
    checkdocs=:exported,
    authors="Davide Ferracin <davide.ferracin@protonmail.com> and contributors",
    pages=["Overview" => "index.md", "Reference" => "reference.md"],
    pagesonly=true,
    linkcheck=true,
    plugins=[bib],
    format=Documenter.HTML(;
        mathengine=Documenter.MathJax(
            Dict(
                :TeX => Dict(
                    :Macros => Dict(
                        :ket => [raw"\lvert #1 \rangle", 1],
                        :bra => [raw"\langle #1 \rvert", 1],
                        :conj => [raw"\overline{#1}", 1],
                        :norm => [raw"\lVert #1\rVert", 1],
                        :abs => [raw"\lvert #1\rvert", 1],
                        :adj => [raw"#1^\dagger", 1],
                        :real => [raw"\operatorname{Re}"],
                        :imag => [raw"\operatorname{Im}"],
                        :sgn => [raw"\operatorname{sgn}"],
                        :N => [raw"\mathbb{N}"],
                        :C => [raw"\mathbb{C}"],
                        :R => [raw"\mathbb{R}"],
                        :dd => [raw"\mathrm{d}"],
                        :env => [raw"_{\mathrm{E}}"],
                        :sys => [raw"_{\mathrm{S}}"],
                        :inter => [raw"_{\mathrm{I}}"],
                        :Max => [raw"_{\mathrm{max}}"],
                    ),
                ),
            ),
        ),
    ),
)

# Automatically deploy documentation to gh-pages.
deploydocs(; repo="github.com/phaerrax/TEDOPA.jl.git")
