from email import header
import pandas as pd
import numpy as np
import os, sys

os.getcwd()
os.chdir("raw")

# ===== ===== ===== ===== ===== Define Functions ===== ===== ===== ===== =====
def format_dt(dt, remove_dup=True):
    dt = dt.copy()

    # rename
    dt.columns = ["id", "chr", "bp", "cM"]

    # remove duplication
    if remove_dup:
        dt["chr"] = dt["chr"].astype(str)
        dt["cM"] = dt["cM"].astype(str)
        dt.loc[:, "tags"] = pd.concat([dt["chr"] + "_" + dt["cM"]])
        dt = dt.loc[~dt["tags"].duplicated()]
        dt = dt.iloc[:, :4]

    # re-type
    dt.loc[:, "chr"] = dt["chr"].astype(int)
    dt.loc[:, "bp"]  = dt["bp"].astype(int)
    dt.loc[:, "cM"]  = dt["cM"].astype(float)

    # sorting
    dt = dt.sort_values(by=["chr", "cM"])
    return dt

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
format_dt(dts).to_csv("../genome_cattle.csv", index=False)

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
format_dt(dts_snp).to_csv("../genome_pig.csv", index=False)

# ===== ===== ===== ===== ===== ===== Maize ===== ===== ===== ===== ===== =====
# Download files for SNP markers that have been genetically mapped
# by at least on map held by MaizeGDB.
# https: // www.maizegdb.org

# Files:
#   MaizeGDB_genetic_SNPs.xref 
#      - SNP ids, position in B73 v3 and v4, associated gene models, 
#        position in 5 major maps (Cornfed Dent Composite, 
#        Cornfed Flint Composite, IBM MaizeSNP50, LHRF Gnp2004, NAM)

#   MaizeGDB_genetic_SNPs.fa
#     - Sequence for SNP markers that are genetically mapped at
#       MaizeGDB.

dt = pd.read_csv("ref_maize.txt", sep="\t", header=None)
# https://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_maize_snp50.pdf
# 55k snp  chip
# b73 v4
# Portwood, J.L., II, Woodhouse, M.R., Cannon, E.K., Gardiner, J.M., Harper, L.C., Schaeffer, M.L., Walsh, J.R., Sen, T.Z., Cho, K.T., Schott, D.A., et al. (2019). MaizeGDB 2018: the maize multi-genome genetics and genomics database. Nucleic Acids Research 47, D1146–D1154.

# select IBM maizeSNP50
dt_s = dt.loc[dt.iloc[:, 10].notna(), :].iloc[:, [0, 2, 3, 11]].dropna() 

# remove string in chromosome column
dt_s.iloc[:, 1] = dt_s.iloc[:, 1].str.replace("Chr", "")

# Export
format_dt(dt_s).to_csv("../genome_maize.csv", index=False)

# ===== ===== ===== ===== ===== ===== Rice ===== ===== ===== ===== ===== =====
# https://rgp.dna.affrc.go.jp/E/publicdata/geneticmap2000/index.html
# Kurata, N., and Yamazaki, Y. (2006). Oryzabase. An integrated biological and genome information database for rice. Plant Physiol 140, 12–17.
dt = pd.read_csv("ref_rice.txt", sep="\t")

dts = dt.iloc[:, [5, 0, 2, 2]]
dts.iloc[:, 2] = 0
format_dt(dts).to_csv("../genome_rice.csv", index=False)


# ===== ===== ===== ===== ===== ===== Maize 3k ===== ===== ===== ===== ===== =====
gd = pd.read_csv("raw_maize_snp.txt", sep="\t")
gm = pd.read_csv("raw_maize_map.txt", sep="\t")

# genotypes
gd = gd.iloc[:, 1:]

# map
gm.columns = ["id", "chr", "bp"]
# re-type
gm.loc[:, "chr"] = dt["chr"].astype(int)
gm.loc[:, "bp"]  = dt["bp"].astype(int)

# export
gd.to_csv("../demo_maize_snp.csv", index=False, header=False)
gm.to_csv("../demo_maize_map.csv", index=False, header=True)

# ===== ===== ===== ===== ===== ===== Verify ===== ===== ===== ===== ===== =====
os.chdir("..")

pd.read_csv("genome_cattle.csv") # 2805
pd.read_csv("genome_pig.csv")    # 2147
pd.read_csv("genome_maize.csv")   # 3746
pd.read_csv("genome_rice.csv")    # 1174

pd.read_csv("demo_maize_snp.csv", header=None) # 281 x 3093
pd.read_csv("demo_maize_map.csv") # 3093





