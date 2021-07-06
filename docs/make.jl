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
        "Basic Usages" => Any[
            "basic/build_genome.md",
            "basic/build_phenome.md",
            "basic/founder.md",
            "basic/mate.md",
            "basic/select.md",
            "basic/breed.md",
        ],
        "Case Studies" => Any[
            "case/simple.md",
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