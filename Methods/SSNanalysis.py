import os
import pandas as pd
from scipy.stats import pearsonr, norm
from collections import defaultdict
from tqdm import tqdm
import shutil
import time
import sys
sys.path.append(r"/proj/c.zihao/work1/1NT/function/")
from SSNfunction import SSNcal


# ====== Configuration ======
#Set dataset name here
dataloop = pd.DataFrame({
    "disease": ["BLCA", "BRCA", "CRC", "ESCA", "HNSC",
                "KIRC", "LIHC", "LUAD", "LUSC", "PRAD", "STAD"],
    "organ": ["Bladder", "Breast", "Colon", "Esophagus", "Salivary Gland",
              "Kidney", "Liver", "Lung", "Lung", "Prostate", "Stomach"]
})


researchAim = '2String9'

for _, row in dataloop.iterrows():
    codestart_time = time.time()

    disease_name = row["disease"]
    organ = row["organ"]
    print(disease_name)

    base_dir = f"/proj/c.zihao/work1/1NT/{researchAim}/{disease_name}"
    save_path = f"{base_dir}/SSN"
    exprSetFile = f"{base_dir}/data/exprSet_filtered.csv"
    NormalFile = "/proj/c.zihao/work1/0ref/GTEx/combined_expr_df.csv"
    link_file = f'/proj/c.zihao/work1/1NT/{researchAim}/links.csv'

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    SSNcal(exprSetFile, NormalFile, link_file, organ, save_path)

    os.chdir(save_path)
    codeend_time = time.time()
    elapsed_time = codeend_time - codestart_time
    print(f"Execution time: {elapsed_time:.4f} seconds")
    with open("runtime_log.txt", "w", encoding="utf-8") as f:
        f.write(f"Execution time: {elapsed_time:.4f} seconds\n")

    print("ALL samples for SSN complete")
    print("=======================================================")