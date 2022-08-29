function DATA(filename::String = ""; header::Bool = true)
    if filename in ["genotypes", "haplotypes", "pedigree", "maize_snp"]
        header = false
    end

    return CSV.read(PATH(filename), DataFrame, header = header)
end


function PATH(filename::String = "")
    root = dirname(dirname(pathof(XSim)))

    if filename == "genotypes"
        return joinpath(root, "data", "demo_genotypes.csv")

    elseif filename == "haplotypes"
        return joinpath(root, "data", "demo_haplotypes.csv")

    elseif filename == "map"
        return joinpath(root, "data", "demo_map.csv")

    elseif filename == "pedigree"
        return joinpath(root, "data", "demo_pedigree.csv")

    # maize data
    elseif filename == "maize_snp"
        return joinpath(root, "data", "demo_maize_snp.csv")

    elseif filename == "maize_map"
        return joinpath(root, "data", "demo_maize_map.csv")

    # reference genome
    elseif filename == "genome_pig"
        return joinpath(root, "data", "genome_pig.csv")

    elseif filename == "genome_cattle"
        return joinpath(root, "data", "genome_cattle.csv")

    elseif filename == "genome_maize"
        return joinpath(root, "data", "genome_maize.csv")

    elseif filename == "genome_rice"
        return joinpath(root, "data", "genome_rice.csv")

    else
        LOG("The available options for DATA() or PATH() are:
        General dataset ['genotypes','haplotypes', 'map', 'pedigree']
        Reference genome ['genome_pig', 'genome_cattle', 'genome_maize', 'genome_rice']
        Maize data ['maize_snp', 'maize_map']", "warn")
        return nothing
    end
end

