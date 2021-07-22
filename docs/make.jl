using Documenter
include("../src/XSim.jl")
using .XSim

makedocs(
    sitename = "XSim.jl",
    doctest  = false,
    clean    = false,
    format   = Documenter.HTML(),
    modules  = [XSim],
    pages = [
        "Home"   => "index.md",
        "Demo: Step by Step"=> "demo.md",
        "Library" => Any[
            "lib/build_genome.md",
            "lib/build_phenome.md",
            "lib/cohort.md",
            "lib/mate.md",
            "lib/select.md",
            "basic/breed.md",
        ],
        "Case Studies" => Any[
            "case/crossbreed.md",
            "case/NAM.md",
        ],
        "Genetic Evaluation" => Any[
            "ge/gblup_pblup.md",
            "ge/multi_trait.md",
        ],
        "For Developers" => Any[
            "dev/custom_data.md",
            "dev/expand_genome.md",
        ],
    ]

)

# deploydocs(
#     repo="github.com/reworkhow/XSim.jl.git",
# )
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

# deploydocs(
#     repo = "github.com/poissonfish/XSim.jl.git",
#     target = "build",
# )