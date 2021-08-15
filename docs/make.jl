using Documenter
include("../src/XSim.jl")
using .XSim

makedocs(
    sitename = "XSim.jl",
    doctest  = false,
    clean    = false,
    format   = Documenter.HTML(
                    analytics="G-RW1CQJ0L6K",
                    collapselevel=4),
    modules  = [XSim],
    pages = [
        "Home"   => "index.md",
        "Demo: Step by Step"=> "demo.md",
        "Core Functions" => Any[
            "core/build_genome.md",
            "core/build_phenome.md",
            "core/cohort.md",
            "core/mate.md",
            "core/select.md",
            "core/GE.md",
            "core/breed.md",
        ],
        "Case Studies" => Any[
            "case/crossbreed.md",
            "case/NAM.md",
        ],
        "Library " => "lib.md",
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
#     repo = "github.com/reworkhow/XSim.jl.git",
#     target = "build",
#     push_preview = true,
# )
