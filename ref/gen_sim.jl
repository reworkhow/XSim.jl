using DelimitedFiles,CSV,Random,DataFrames,LinearAlgebra,SparseArrays,Statistics,StatsBase
Pkg.add("XSim")
using XSim
map = CSV.read("/Users/lijinghui/Desktop/Genetics/Group_GWAS/pig_data/map_rm_minor.csv", DataFrame)

# number of loci for each chromosome
numLoci = [0]
for i in 1:18
    chr = findall(map.chrs .== i)
    push!(numLoci, length(chr))
end
numLoci = numLoci[2:end]

# length of each chromosome (Vingborg et al., 2009)
chrLength = [151.4, 150.7, 100.5, 137.6, 53.4, 148, 64.2, 107.4, 116.6, 139.5, 15.7, 69.8, 128, 96.3, 100.3, 38.1, 60.2, 34.1]

# number of chromosome; mutation rate
numChr,rate_mutation = 18,0.0

# gene frequency, 0.5 for all
geneFreq = [[0.5]]
for i in 1:18
    append!(geneFreq, [ones(numLoci[i])*0.5])
end
geneFreq = geneFreq[2:end]

for replicate in 1:10
    Random.seed!(100 + replicate)
    # QTL position sumulation
    qtl_pos = rand(1:size(map)[1], 100)
    sort!(qtl_pos)
    qtl_chr = map.chrs[qtl_pos]
    n_qtl = []
    for i in 1:18
        chr = findall(qtl_chr .== i)
        push!(n_qtl, length(chr))
    end

    # QTL index (corresponds to loci index)
    lociIndex0 = [0]
    for i in 1:18
        append!(lociIndex0, Vector(1:numLoci[i]))
    end
    lociIndex = lociIndex0[2:end]
    qtl_pos0 = lociIndex[qtl_pos]

    qtlIndex = [[1]]
    k = 1
    for i in 1:18
        append!(qtlIndex, [qtl_pos0[k:sum(n_qtl[1:i])]])
        k = k + n_qtl[i]
    end
    qtlIndex = qtlIndex[2:end]

    # map position of loci in bp
    mapPos0 = map.pos
    mapPos = [[0.0]]
    k = 1
    for i in 1:18
        append!(mapPos, [mapPos0[k:sum(numLoci[1:i])]])
        k = k + numLoci[i]
    end
    mapPos = mapPos[2:end]

    # map position of loci in cM
    mapPos_cM = [[0.0]]
    for i in 1:18
        append!(mapPos_cM, [mapPos[i] .* (chrLength[i] - 0.05) ./ mapPos[i][end]]) # -0.05 so that last locus is a bit shorter than chorm length
    end


    # QTL effects simulation
    qtlEffects = [[0.0]]
    for i in 1:18
        append!(qtlEffects, [randn(n_qtl[i])])
    end
    qtlEffects = qtlEffects[2:end]
    QTL_dat = DataFrame(QTL_pos = qtl_pos, QTL_eff = vcat(qtlEffects...))
    CSV.write("/Users/lijinghui/Desktop/Server_result/Group_GWAS/simple_sim/"*string(replicate)*"/QTL_dat.csv", QTL_dat)

    # Build founders based on known haplotypes
    build_genome(numChr,chrLength,numLoci,geneFreq,mapPos_cM,qtlIndex,qtlEffects)
    sireSizeFounder = 30
    damSizeFounder = 20
    sires = sampleFounders(sireSizeFounder,"/Users/lijinghui/Desktop/Genetics/Group_GWAS/pig_data/male_founder_hap_wo_missing.txt")
    dams = sampleFounders(damSizeFounder,"/Users/lijinghui/Desktop/Genetics/Group_GWAS/pig_data/female_founder_hap_wo_missing.txt")

    sire_geno = float(getOurGenotypes(sires))
    dam_geno = float(getOurGenotypes(dams))

    nSires,nDams = 30,20
    popSize,ngen = 500,2
    varRes = 0.0

    sire2,dam2,gen2=sampleSel(popSize, nSires, nDams, ngen, sires, dams, varRes);

    sireID = []
    damID = []
    for i in 1:250
        push!(sireID, sire2.animalCohort[i].sireID)
        push!(damID, sire2.animalCohort[i].damID)
    end

    sire_selected = sort(unique(sireID))
    dam_selected = sort(unique(damID))

    qtl_eff = []
    for i in 1:18
        append!(qtl_eff, qtlEffects[i])
    end
    qtl_whole = zeros(size(sire_geno)[2])
    qtl_whole[qtl_pos] = qtl_eff

    sire2_geno = float(getOurGenotypes(sire2))
    dam2_geno = float(getOurGenotypes(dam2))
    geno_gen2 = vcat(sire2_geno, dam2_geno)
    writedlm("/Users/lijinghui/Desktop/Server_result/Group_GWAS/simple_sim/"*string(replicate)*"/genotype.txt", geno_gen2)

end
