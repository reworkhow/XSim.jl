import pandas as pd
import numpy as np
import os, sys

os.getcwd()

# ===== ===== ===== ===== ===== ===== Cattle ===== ===== ===== ===== ===== =====
dt = pd.read_csv("ref_cattle.csv")

# Filter dataset
dts = dt.loc[
        # select 1 ~ 29 chromosome
        (dt["Chr"].isin([str(i + 1) for i in range(29)])) &
        # select SNP variants
        (dt["Type"] == "SNP")].\
        iloc[:, [3, 0, 8, 7]].\
        dropna()

# Format dataset
dts.columns = ["id", "chr", "bp", "cM"]
dts.chr = dts.chr.astype(int)
dts.bp  = dts.bp.astype(int)

# Export
dts.to_csv("genome_cattle.csv", index=False)

# ===== ===== ===== ===== ===== ===== Pig ===== ===== ===== ===== ===== =====
# SNP dataset
dt_snp = pd.read_csv("ref_pig_60k.csv", skiprows=9, sep="\t")
dts_snp = dt_snp.loc[
            # select 1 ~ 18 chromosome
            (dt_snp["Chromosome"].isin([str(i + 1) for i in range(18)]))].\
            iloc[:, [0, 1, 2]]
dts_snp = dts_snp.reset_index().iloc[:, 1:]
dts_snp.columns = ["id", "chr", "bp"]
dts_snp.chr = dts_snp.chr.astype(int)

# Linkage dataset
dt_lkg = pd.read_csv("ref_pig_linkage.txt", skiprows=4, sep="\t")
dts_lkg = dt_lkg.iloc[:, [0, 1, 3, 5, 6]]
dts_lkg.columns = ["chr", "bin", "bp_start", "bp_end", "rec_rate"]
dts_lkg.loc[:, "rec_rate"] = dts_lkg["rec_rate"].str.replace(',', '.').astype(float)
dts_lkg.loc[:, "cM"] =  dts_lkg.loc[:, ["chr", "bin", "rec_rate"]].\
                            groupby(["chr", "bin"]).sum().\
                            groupby(level=0).cumsum().\
                            reset_index().\
                            loc[:, "rec_rate"].values

# Merge cM to SNP dataset
np_cM = np.zeros(len(dts_snp))
for idx, row in dts_snp.iterrows():
    cM = dts_lkg.loc[(dts_lkg.bp_start <= row.bp) &
                     (dts_lkg.bp_end   >= row.bp) &
                     (dts_lkg.chr      == row.chr), :].cM.values
    if len(cM) != 0:
        np_cM[idx] = round(cM[0], 3)

dts_snp.loc[:, "cM"] = np_cM
dts_snp = dts_snp.loc[dts_snp.cM > 0]


# Export
dts_snp.sort_values(["chr", "cM"]).to_csv("genome_pig.csv", index=False)


# ===== ===== ===== ===== ===== ===== Verify ===== ===== ===== ===== ===== =====
pd.read_csv("genome_cattle.csv") # 6,231
pd.read_csv("genome_pig.csv")    # 45,292

