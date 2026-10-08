using Documenter
using NetworkCore
using Graphs   # for `Graphs.is_bipartite`, which NetworkCore.jl extends but does not export

DocMeta.setdocmeta!(NetworkCore, :DocTestSetup, :(using NetworkCore); recursive=true)

makedocs(
    sitename = "NetworkCore.jl",
    modules = [NetworkCore],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "DOCS_PRETTY_URLS", get(ENV, "CI", "false")) == "true",
        canonical = "https://statistical-network-analysis-with-Julia.github.io/NetworkCore.jl/dev/",
        assets = ["assets/favicon.ico"],
        footer = "[Ecosystem home](/) · [Packages](/packages/) · [Get started](/getting-started/) · [Capabilities](/capabilities/) — Built with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
        edit_link = "main",
    ),
    repo = Documenter.Remotes.GitHub("Statistical-network-analysis-with-Julia", "NetworkCore.jl"),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Creating Networks" => "guide/construction.md",
            "Attributes" => "guide/attributes.md",
            "Missing Data" => "guide/missing_data.md",
            "Conversion" => "guide/conversion.md",
            "Conversion Invariants" => "guide/conversion_invariants.md",
            "I/O" => "guide/io.md",
            "Datasets" => "guide/datasets.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Graph Interface" => "api/graph_interface.md",
            "Attributes" => "api/attributes.md",
            "Conversion & I/O" => "api/conversion.md",
            "Ecosystem Contracts" => "api/contracts.md",
            "Result Metadata" => "api/metadata.md",
            "Shared Inference" => "api/inference.md",
        ],
    ],
    # STRICT. Undefined bindings, bad cross-references, duplicate docs and
    # malformed markdown are build ERRORS, so they cannot silently accumulate
    # again (a docs build that passes while warning is one that will rot).
    #
    # `checkdocs = :exports` is the one deliberate exclusion: every *exported*
    # name must be documented, but internal machinery need not be. That covers
    # the materialized/internal types a user never names — `UndirectedEdgeIter`
    # (the undirected `edges()` iterator), the `Base.copy`/`SimpleGraph`/
    # `SimpleDiGraph`/`Graphs.is_bipartite` method extensions, and the inner
    # constructors — which would otherwise demand filler docstrings.
    warnonly = false,
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/statistical-network-analysis-with-Julia/NetworkCore.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev", # Development alias; change to "v^" when adopting release-based stable docs.
        "dev" => "dev",
    ],
    push_preview = false, # Pull requests build docs; main/tags publish through Pages.
)
