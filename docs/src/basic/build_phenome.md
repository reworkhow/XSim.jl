# Build Phenome

## Data

```
# map.csv
id,chr,bp,cM,MAF,eff_1,eff_2
snp 1,1,1818249,50.8,0.5,1.5,2.8
snp 2,1,6557697,80.3,0.5,0.0,0.0
snp_3,2,2298800,39.2,0.5,0.0,0.0
snp 4,2,5015698,66.3,0.5,0.0,0.0
```

## Assign by file

```julia;
build_phenome("map.csv", h2 = 0.3)
```

## Assign by n_qtl

```julia
build_phenome(n_qtl = 100,
              Vg    = [1 .5
                      .5  1],
              h2    = [0.3, 0.8])
```