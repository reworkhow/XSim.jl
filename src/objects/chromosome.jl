abstract type AbstractAnimal end

mutable struct Chromosome
    index::Int64
    haplotype::Array{AlleleIndexType,1}
    ori::Array{Int64,1}
    pos::Array{Float64,1}
    mut::Array{Float64,1}

    # Founder's chromosome
    function Chromosome(i_chr::Int64,
        ori::Int64,
        haplotypes::Array{AlleleIndexType,1} = [0])

        n_loci = GLOBAL("n_loci", chromosome = i_chr)

        if haplotypes != [0]
            # Alleles provided by files
            chromosome = new(i_chr,
                haplotypes,
                [ori],
                [0.0],
                Array{Float64}(undef, 0))

        else
            # Alleles provided by Bernoulli sampler
            chromosome = new(i_chr,
                Array{AlleleIndexType}(undef, n_loci),
                [ori],
                [0.0],
                Array{Float64}(undef, 0))
            set_alleles!(chromosome, n_loci)
        end

        return chromosome
    end

    function set_alleles!(chromosome::Chromosome,
        n_loci::Int)

        i_chr = chromosome.index
        i_maf = GLOBAL("maf", chromosome = i_chr)
        for j_locus in 1:n_loci
            p = Bernoulli(1 - i_maf[j_locus])
            chromosome.haplotype[j_locus] = convert(AlleleIndexType, rand(p))
        end
    end


    # Progeny's chromosome
    function Chromosome(i_chr::Int64,
        parent::AbstractAnimal)

        # n_loci = GLOBAL("n_loci", chromosome = i_chr)
        chromosome = new(i_chr,
            Array{AlleleIndexType}(undef, 0),
            Array{Int64}(undef, 0),
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0))
        sample_genome!(chromosome, parent)

        return chromosome
    end
end

Base.show(io::IO, chromosome::Chromosome) = print(chromosome)
function print(chromosome::Chromosome)
    println("index     : ", chromosome.index)
    println("haplotype : ", chromosome.haplotype)
    println("origin    : ", chromosome.ori)
    println("position  : ", chromosome.pos)
    println("mutation  : ", chromosome.mut)
end

function set_haplotypes!(genome::Array{Chromosome,1})

    for i_chr in 1:GLOBAL("n_chr")
        numLoci = GLOBAL("n_loci", chromosome = i_chr)
        resize!(genome[i_chr].haplotype, numLoci)

        numOri = length(genome[i_chr].ori)
        push!(genome[i_chr].pos, GLOBAL("length_chr", chromosome = i_chr))

        i_locus = 1
        position = GLOBAL("cM", chromosome = i_chr, locus = i_locus) / 100

        lociPerM = round(Int64, numLoci / genome[i_chr].pos[numOri+1])
        segLen = 0
        prevSegLen = 0
        endLocus = 0

        for segment in 1:numOri
            whichFounder = ceil(Integer, genome[i_chr].ori[segment] / 2)
            genomePatorMatInThisFounder = (genome[i_chr].ori[segment] % 2 == 0) ?
                                          GLOBAL("founders")[whichFounder].genome_dam[i_chr] :
                                          GLOBAL("founders")[whichFounder].genome_sire[i_chr]

            startPos = genome[i_chr].pos[segment]
            endPos = genome[i_chr].pos[segment+1]
            prevSegLen += segLen
            segLen = 0
            if segment < numOri
                numLociUntilGuessedPos = round(Int64, endPos * lociPerM)
                numLociUntilGuessedPos = maximum([1, numLociUntilGuessedPos])
                if numLociUntilGuessedPos > numLoci
                    numLociUntilGuessedPos = numLoci
                end

                guessedPos = GLOBAL("cM",
                    chromosome = i_chr,
                    locus = numLociUntilGuessedPos) / 100
                if guessedPos > endPos
                    ul = numLociUntilGuessedPos
                    ll = i_locus
                    if ll > endPos
                        ll -= 1
                    end
                else
                    ll = numLociUntilGuessedPos
                    ul = numLoci
                end

                iter = 0
                while ul - ll > 1
                    prevNumLociUntilGuessedPos = numLociUntilGuessedPos
                    iter += 1
                    numLociUntilGuessedPos = (ul - ll) / 2
                    numLociUntilGuessedPos = ll + round(Int64, numLociUntilGuessedPos)

                    guessedPos = GLOBAL("cM",
                        chromosome = i_chr,
                        locus = numLociUntilGuessedPos) / 100
                    if prevNumLociUntilGuessedPos == numLociUntilGuessedPos
                        if guessedPos == ll
                            pos_ll = GLOBAL("cM", chromosome = i_chr, locus = ll + 1) / 100
                            if pos_ll < endPos
                                ll += 1
                            else
                                ul -= 1
                            end
                        else
                            pos_ul = GLOBAL("cM", chromosome = i_chr, locus = ul - 1) / 100
                            if pos_ul > endPos
                                ul -= 1
                            else
                                ll += 1
                            end
                        end
                    elseif guessedPos > endPos
                        ul = numLociUntilGuessedPos
                    else
                        ll = numLociUntilGuessedPos
                    end
                end
                endLocus = ll
                segLen = ll - prevSegLen
                if i_locus < numLoci
                    i_locus = ll + 1
                end
            elseif i_locus <= numLoci
                segLen = numLoci - endLocus
                endLocus = numLoci
            else
                segLen = 0
            end
            if segLen > 0
                genome[i_chr].haplotype[(endLocus-segLen+1):endLocus] = genomePatorMatInThisFounder.haplotype[(endLocus-segLen+1):endLocus]
            end
        end

        for j in 1:length(genome[i_chr].mut)
            whichlocus = findfirst((GLOBAL("cM", chromosome = i_chr) ./ 100) .== genome[i_chr].mut[j])
            genome[i_chr].haplotype[whichlocus] = 1 - genome[i_chr].haplotype[whichlocus]
        end

        pop!(genome[i_chr].pos) #remove temporary added chrLength
    end
end