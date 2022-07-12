# XSim.jl

Documentation for XSim.jl

!!! tip "New Users"
    New users are strongly encouraged to read the page [Demo: Step by Step](@ref) first.

## Features

- Efficient CPOS algorithm
- Founders characterized by real genome sequence data
- Complicated pedigree structures among descendants
- Strong extensibilty
- Modularism design
- Straightforward interface

## Installation
```jldoctest
julia> using Pkg
julia> Pkg.add("XSim")
```
or the beta version
```jldoctest
julia> Pkg.add(PackageSpec(name="XSim", rev="master"))
```

## Outline
##### Home
```@contents
Pages = ["index.md"]
```
##### Demo
```@contents
Pages = ["demo.md"]
```
##### Core Functions
```@contents
Pages = ["core/build_genome.md"]
```
```@contents
Pages = ["core/build_phenome.md"]
```
```@contents
Pages = ["core/cohort.md"]
```
```@contents
Pages = ["core/mate.md"]
```
```@contents
Pages = ["core/select.md"]
```
```@contents
Pages = ["core/GE.md"]
```
```@contents
Pages = ["core/breed.md"]
```
##### Case Studies
```@contents
Pages = ["case/crossbreed.md"]
```
```@contents
Pages = ["case/NAM.md"]
```
##### Library
```@contents
Pages = ["lib.md"]
```

## Cite XSimV2
**BibTex**
```BibTeX
@article{chen_xsim_2022,
 title = {XSim version 2: simulation of modern breeding programs},
 volume = {12},
 issn = {2160-1836},
 url = {<https://doi.org/10.1093/g3journal/jkac032>},
 doi = {10.1093/g3journal/jkac032},
 number = {4},
 urldate = {2022-05-26},
 journal = {G3 Genes{\textbar}Genomes{\textbar}Genetics},
 author = {Chen, Chunpeng James and Garrick, Dorian and Fernando, Rohan and Karaman, Emre and Stricker, Chris and Keehan, Michael and Cheng, Hao},
 month = apr,
 year = {2022},
}
```
**Bibliography**
> Chen, C.J., D. Garrick, R. Fernando, E. Karaman, C. Stricker, M. Keehan, and H. Cheng. 2022. XSim version 2: simulation of modern breeding programs. G3 Genes|Genomes|Genetics 12:jkac032. doi:10.1093/g3journal/jkac032.
