function sampleOnePosOri(genome::Array{Chromosome,1}, parent::Animal)
    numberChromosomePair = get_num_chrom(common.G)

    for i in 1:numberChromosomePair

        genome[i] = Chromosome(i,
                               Array{Int64  }(undef, 0),
                               Array{Int64  }(undef, 1),
                               Array{Float64}(undef, 1),
                               Array{Float64}(undef, 1))

        currentChrom = (rand(Bernoulli(0.5)) == 1) ? parent.genome_sire[i] : parent.genome_dam[i]

        chrLength = common.G.chr[i].chrLength

        binomialN = convert(Int64, ceil(chrLength * 3 + 1))
        numCrossover = rand(Binomial(binomialN, chrLength / binomialN))
        rec = [0.0]

        for irec in 1:numCrossover
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

            currentChrom = (currentChrom == parent.genome_dam[i]) ? parent.genome_sire[i] : parent.genome_dam[i]

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

        resize!(genome[i].pos, numTemp)
        resize!(genome[i].ori, numTemp)

        #add new mutation (no merge required as in pos/ori below)
        resize!(genome[i].mut, numTempMut)

        for muti in 1:numTempMut
            genome[i].mut[muti] = tempMut[muti]
        end

        mutation_rate = common.G.mutRate
        numLoci       = common.G.chr[i].numLoci
        nmut          = rand(Binomial(numLoci, mutation_rate))
        if nmut != 0
            muts = sample(common.G.chr[i].mapPos, nmut)
            genome[i].mut = vcat(genome[i].mut, muts)
            sort!(genome[i].mut)
        end

        ##merging of nearby same ori;might be more effcient

        genome[i].pos[1] = tempPos[1]
        genome[i].ori[1] = tempOri[1]

        this = 1
        for m in 2 : numTemp
            if(tempOri[m] != tempOri[m - 1])
                this = this + 1
                genome[i].pos[this] = tempPos[m]
                genome[i].ori[this] = tempOri[m]
            end
        end

        resize!(genome[i].pos,this) #resize after merging
        resize!(genome[i].ori,this)
    end
end
