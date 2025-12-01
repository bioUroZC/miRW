import os
import time
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import shutil
import sys
sys.path.append(r"/proj/c.zihao/work1/1NT/function/")
from SWEETfunction import SWEETcal

# ------------------------- Set input paths -------------------------

#Set dataset name here
available_datasets = ["BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KIRC",
                       "LIHC", "LUAD", "LUSC", "PRAD", "STAD"]

researchAim = '2String9'

for disease_name in available_datasets:
    # === Setup Working Directory and Load Global Data ===
    print(disease_name)
    base_path = f"/proj/c.zihao/work1/1NT/{researchAim}/{disease_name}"
    exprSetFile = f"{base_path}/data/exprSet_filtered.csv"
    save_path = f"{base_path}/SWEET"
    link_file = f'/proj/c.zihao/work1/1NT/{researchAim}/links.csv'

    codestart_time = time.time()
    print(codestart_time)

    if os.path.exists(save_path):
        shutil.rmtree(save_path)
    os.makedirs(save_path)

    SWEETcal(exprSetFile, link_file, save_path)

    os.chdir(save_path)
    codeend_time = time.time()
    elapsed_time = codeend_time - codestart_time
    print(f"Execution time: {elapsed_time:.4f} seconds")
    with open("runtime_log.txt", "w", encoding="utf-8") as f:
        f.write(f"Execution time: {elapsed_time:.4f} seconds\n")


    print("[INFO] All sameple for SWEET saved.")
    print("=======================================================")
