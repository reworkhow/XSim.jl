using Documenter
using XSim

makedocs(
    sitename = "XSim.jl",
    doctest  = false,
    clean    = false,
    format   = Documenter.HTML(),
    modules  = [XSim],
    pages = [
        "Home"   => "index.md",
        "Getting Started"=> "demo.md",
        "Library" => Any[
            "lib/build_genome.md",
            # "lib/build_phenome.md",
            # "lib/founder.md",
            # "lib/mate.md",
            # "lib/select.md",
            # "lib/breed.md",
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