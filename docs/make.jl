using Documenter
using XSim

makedocs(
    sitename = "XSim",
    doctest  = false,
    clean    = true,
    format   = Documenter.HTML(),
    modules  = [XSim],
    pages = Any[
        "Home"   => "index.md",
        "Basic Usages" => Any[
            "basic/build_genome.md",
            "basic/build_phenome.md",
            "basic/founder.md",
            "basic/mate.md",
            "basic/select.md",
            "basic/breed.md",
        ],
        "Case Studies" => Any[
            "case/crossbreed.md",
            "case/NAM.md",
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