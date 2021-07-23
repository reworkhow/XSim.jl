function genetic_evaluation(cohort          ::Cohort,
                            phenotypes     ::DataFrame=DataFrame();
                            model_equation ::String="",
                            covariate      ::String="",
                            random_iid     ::String="",
                            random_str     ::String="",
                            methods        ::String="GBLUP",
                            idx_missing    ::Any=[],
                            return_out     ::Bool=false,
                            args...)

    # 1. Acquire phenotypes
    if nrow(phenotypes) == 0
        phenotypes = get_phenotypes(cohort, "JWAS"; args...)
    end

    # 1.1 If missing is provided
    if length(idx_missing) != 0
        XSim.allowmissing!(phenotypes);
        phenotypes[idx_missing, 2:end] .= missing
    end

    # 2. Build Model Equations
    if model_equation == ""
        traits = names(phenotypes)[2:end]
        array_eq = ["$(trait) = intercept" for trait in traits]
        model_equation = join(array_eq, "\n")
    end

    # 3. Build model
    model = JWAS.build_model(model_equation);

    # 4: Set Factors or Covariates
    if covariate != ""
        JWAS.set_covariate(model, covariate);
    end

    # 5: Set Random or Fixed Effects
    if random_iid != ""
        JWAS.set_random(model, random_iid);
    end
    if random_str != ""
        pedigree = get_pedigree(cohort,   "JWAS")
        JWAS.set_random(model, random_str, pedigree);
    end

    # 6. If GBLUP, add genotypes
    if methods == "GBLUP"
        # Add genotype for GBLUP
        genotypes= get_genotypes(cohort) # 0 1 2
        JWAS.add_genotypes(model, float.(genotypes))
    end

    # 7. Run MCMC
    out = JWAS.runMCMC(model, phenotypes, methods=methods);

    # Remove outputs
    try
        rm("results", recursive=true)
        rm("IDs_for_individuals_with_genotypes.txt")
        rm("IDs_for_individuals_with_pedigree.txt")
        rm("IDs_for_individuals_with_phenotypes.txt")
    catch e
        nothing
    end

    # Outputs
    if return_out
        # return complete JWAS report
        return out
    else
        # return EBV only
        return hcat([out["EBV_$x"][:, "EBV"] for x in traits]...)
    end

    # Note
    # jwas_ped      = get_pedigree(cohort, "JWAS");
    # jwas_p        = get_phenotypes(cohort, "JWAS");
    # allowmissing!(jwas_p);
    # true_p = jwas_p.y1[1:200]
    # jwas_p.y1[1:200] .= missing;

    # genotypes     = get_genotypes(cohort, "JWAS");
    # model         = XSim.JWAS.build_model(model_equation);
    # out = XSim.JWAS.runMCMC(model, jwas_p, methods="GBLUP")

    # Not working
    # eq="y1 = intercept + genotypes"
    # genetic_evaluation(cohort, model_equation=eq)
end
