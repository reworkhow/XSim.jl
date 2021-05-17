function sample_genome!(chromosome::Chromosome, parent::Animal)

    i = chromosome.index

    # randomly select sire or dam genome
    currentChrom = (rand(Bernoulli(0.5)) == 1) ? parent.genome_sire[i] : parent.genome_dam[i]

    chrLength = GLOBAL("length_chr", chromosome=i) #### it's Morgan, not cM!

    n_Binomial = convert(Int64, ceil(chrLength * 3 + 1))
    p_Binomial = convert(Float64, chrLength / n_Binomial)
    numCrossover = rand(Binomial(n_Binomial, p_Binomial))
    rec = [0.0]

    for i in 1:numCrossover
        push!(rec, chrLength * rand(1)[1])
    end

    push!(rec, chrLength)
    sort!(rec)                #rec is like 0.00,0.23,0.45,0.76,1.00

    numTemp = 1
    numTempMut = 0
    for j in 1:(length(rec) - 1)
        for k in 1:length(currentChrom.pos)
            if currentChrom.pos[k] >= rec[j] && currentChrom.pos[k] < rec[j + 1]
                tempPos[numTemp] = currentChrom.pos[k]
                tempOri[numTemp] = currentChrom.ori[k]
                numTemp = numTemp + 1
            elseif currentChrom.pos[k] >= rec[j + 1]
                break
            end
        end
        for k in 1:length(currentChrom.mut)
            if currentChrom.mut[k] >= rec[j] && currentChrom.mut[k] < rec[j + 1]
                numTempMut = numTempMut + 1
                tempMut[numTempMut] = currentChrom.mut[k]
            elseif currentChrom.mut[k] >= rec[j + 1]
                break
            end
        end

        currentChrom = (currentChrom == parent.genome_dam[i]) ?
                            parent.genome_sire[i] :
                            parent.genome_dam[i]

        findRecOri = 0
        m = 1
        lengthChrPos = length(currentChrom.pos)
        while m <= lengthChrPos && currentChrom.pos[m] <= rec[j + 1] #pos[m] cannot be the length of chromosome
        m += 1
        findRecOri += 1
        end
        tempPos[numTemp] = rec[j + 1]
        tempOri[numTemp] = currentChrom.ori[findRecOri]
        numTemp = numTemp + 1 #**
    end

    numTemp = numTemp - 1 #remove the last one got from **
    numTemp = numTemp - 1 #remove the last one which is the length of chromsome

    resize!(chromosome.pos, numTemp)
    resize!(chromosome.ori, numTemp)

    #add new mutation (no merge required as in pos/ori below)
    resize!(chromosome.mut, numTempMut)

    for muti in 1:numTempMut
        chromosome.mut[muti] = tempMut[muti]
    end

    mutation_rate = GLOBAL("rate_mutation")
    numLoci       = GLOBAL("n_loci", chromosome=i)
    nmut          = rand(Binomial(numLoci, mutation_rate))
    if nmut != 0
        muts = sample(GLOBAL("cM", chromosome=i), nmut) ./ 100
        chromosome.mut = vcat(chromosome.mut, muts)
        sort!(chromosome.mut)
    end

    ##merging of nearby same ori;might be more effcient

    chromosome.pos[1] = tempPos[1]
    chromosome.ori[1] = tempOri[1]

    this = 1
    for m in 2 : numTemp
        if(tempOri[m] != tempOri[m - 1])
            this = this + 1
            chromosome.pos[this] = tempPos[m]
            chromosome.ori[this] = tempOri[m]
        end
    end

    resize!(chromosome.pos, this)
    resize!(chromosome.ori, this)
end