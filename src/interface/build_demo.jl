function build_demo()
    n_chr = 10
    n_loci_chr = 100
    n_loci = n_chr * n_loci_chr

    chromosome = [i         for i in 1:n_chr for j in 1:n_loci_chr]
    bp         = [10 * j    for i in 1:n_chr for j in 1:n_loci_chr]
    cM         = [1.5 * j   for i in 1:n_chr for j in 1:n_loci_chr]
    maf        = fill(0.5, n_loci)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf;
                 rate_mutation=rate_mutation,
                 rate_error=rate_error)

    n_qtl = [3, 8]
    vg    = [1.0  0
             0  1.0]
    build_phenome(n_qtl; vg=vg)

end

function build_demo_small()
    println("This is a demo data set")
    n_chr = 2
    n_loci_chr = 5
    n_loci = n_chr * n_loci_chr

    chromosome = [i         for i in 1:n_chr for j in 1:n_loci_chr]
    bp         = [10 * j    for i in 1:n_chr for j in 1:n_loci_chr]
    cM         = [1.5 * j   for i in 1:n_chr for j in 1:n_loci_chr]
    maf        = fill(0.5, n_loci)
    rate_mutation = 0.0
    rate_error    = 0.0
    build_genome(chromosome, bp, cM, maf;
                 rate_mutation=rate_mutation,
                 rate_error=rate_error)

    n_qtl = [2, 4]
    vg    = [1.0  0
             0  1.0]
    build_phenome(n_qtl; vg=vg)

end