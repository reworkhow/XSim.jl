# Genetic Evaluation
Users can select individuals by random, phenotypes, or esti- mated breeding values from genetic evaluations. A genome-enabled analysis package [JWAS](https://github.com/reworkhow/JWAS.jl) (Cheng et al. 2018) has been already incorporated into XSimV2. Multiple methods can be per- formed in XSimV2 for genetic evaluations, including pedigree- based BLUP (Henderson 1984), GBLUP (Habier et al. 2007; Van- Raden 2008), Bayesian Alphabet (Meuwissen et al. 2001; Park and Casella 2008; Kizilkaya et al. 2010; Habier et al. 2011; Erbe et al. 2012; Moser et al. 2015; Gianola and Fernando 2020), and single-step methods (Legarra et al. 2009; Fernando et al. 2014) for both single-trait and multiple-trait analysis (Gianola and Fernando 2020).

```@docs
genetic_evaluation
```

## References
Cheng, H., R. Fernando, and D. Garrick, 2018 JWAS: Julsi- taaimnoptlheemrentation of Whole-genome Analysis Soft- ware. Proceedings of the world congress on genetics applied to livestock production 11: 859.

Erbe, M., B. J. Hayes, L. K. Matukumalli, S. Goswami, P. J. Bow- man, et al., 2012 Improving accuracy of genomic predictions within and between dairy cattle breeds with imputed high- density single nucleotide polymorphism panels. Journal of Dairy Science 95: 4114–4129.

Fernando, R. L., J. C. Dekkers, and D. J. Garrick, 2014 A class of Bayesian methods to combine large numbers of genotyped and non-genotyped animals for whole-genome analyses. Ge- netics Selection Evolution 46: 50.

Gianola, D. and R. L. Fernando, 2020 A Multiple-Trait Bayesian Lasso for Genome-Enabled Analysis and Prediction of Com- plex Traits. Genetics 214: 305–331, Publisher: Genetics Section: Investigations.

Habier, D., R. L. Fernando, and J. C. M. Dekkers, 2007 The Im- pact of Genetic Relationship Information on Genome-Assisted Breeding Values. Genetics 177: 2389–2397, Publisher: Genetics Section: Investigations.

Habier, D., R. L. Fernando, K. Kizilkaya, and D. J. Garrick, 2011 Extension of the bayesian alphabet for genomic selection. BMC Bioinformatics 12: 186.

Henderson, C. R., 1984 Applications of linear models in animal breeding.. Publisher: University of Guelph.

Kizilkaya, K., R. L. Fernando, and D. J. Garrick, 2010 Genomic prediction of simulated multibreed and purebred performance using observed fifty thousand single nucleotide polymor- phism genotypes. Journal of Animal Science 88: 544–551.

Legarra, A., I. Aguilar, and I. Misztal, 2009 A relationship matrix including full pedigree and genomic information. Journal of Dairy Science 92: 4656–4663.

Meuwissen, T. H. E., B. J. Hayes, and M. E. Goddard, 2001 Pre- diction of Total Genetic Value Using Genome-Wide Dense Marker Maps. Genetics 157: 1819–1829, Publisher: Genetics Section: Investigations.

Moser, G., S. H. Lee, B. J. Hayes, M. E. Goddard, N. R. Wray, et al., 2015 Simultaneous Discovery, Estimation and Prediction Analysis of Complex Traits Using a Bayesian Mixture Model. PLOS Genetics 11: e1004969, Publisher: Public Library of Science.

Park, T. and G. Casella, 2008 The Bayesian Lasso. Journal of the American Statistical Association 103: 681–686.

VanRaden, P. M., 2008 Efficient Methods to Compute Genomic
Predictions. Journal of Dairy Science 91: 4414–4423.
