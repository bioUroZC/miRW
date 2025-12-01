#!/bin/bash
#SBATCH -J 1PPP
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH -o 1PPP.out
#SBATCH -e 1PPP.err

start_time=$(date +%s)
start_str=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job started at: $start_str" >> 1PPP.out

cd /proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/PPIXpress123
java -jar PPIXpress.jar -t=3 /proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/human_ppin.sif.gz /proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/out /proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BL_A13J_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BL_A13J_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20N_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20N_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20Q_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20Q_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20R_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20R_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20U_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20U_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20W_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A20W_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A2LA_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A2LA_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A2LB_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_BT_A2LB_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_CU_A0YN_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_CU_A0YN_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_CU_A0YR_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_CU_A0YR_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A3BM_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A3BM_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A3WC_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A3WC_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A6I3_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GC_A6I3_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A2C5_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A2C5_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A3OP_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A3OP_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A3OQ_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_GD_A3OQ_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A3WV_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A3WV_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A54R_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A54R_11A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A5RI_01A.gz \
/proj/c.zihao/work1/1NT/2String9/BLCA/PPIX/files/TCGA_K4_A5RI_11A.gz

end_time=$(date +%s)
end_str=$(date +"%Y-%m-%d %H:%M:%S")
echo "Job finished at: $end_str" >> 1PPP.out

elapsed=$((end_time - start_time))
echo "Total execution time: ${elapsed} seconds" >> 1PPP.out
