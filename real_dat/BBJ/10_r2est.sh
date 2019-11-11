#!/bin/bash
#SBATCH --partition=nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=test
#SBATCH --mem=2G

#SBATCH --array=1-44
#SBATCH --output=/net/mulan/disk2/yasheng/out/SH_test_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/SH_test_%a.err

SHtest=/net/mulan/home/yasheng/Biobank/code/out_bbj/10_r2est.R

bash
let k=0

for method in ldpredinf ldpredpbest ;do
for PHENO in 10;do
for cross in 4 ; do
for ((chr=1;chr<23;chr++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript ${SHtest} --pheno ${PHENO} --method ${method} --cross ${cross} --chr ${chr}
fi
done
done
done
done
