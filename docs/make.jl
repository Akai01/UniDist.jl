using Documenter
using UniDist

makedocs(
    sitename = "UniDist.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Akai01.github.io/UniDist.jl",
    ),
    modules = [UniDist],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Distributions" => [
            "distributions/continuous.md",
            "distributions/discrete.md",
        ],
        "API Reference" => [
            "api/core.md",
            "api/survival.md",
            "api/intervals.md",
        ],
        "Examples" => [
            "examples/basic-usage.md",
            "examples/vectorized.md",
            "examples/statistical-analysis.md",
            "examples/bayesian.md",
        ],
    ]
)

deploydocs(
    repo = "github.com/Akai01/UniDist.jl.git",
    devbranch = "master",
)
