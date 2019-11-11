#!/bin/bash

#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lassosum
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

#SBATCH --array=1-3
#SBATCH --output=/net/mulan/disk2/yasheng/out/lassosum%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/lassosum%a.err

lassosum=/net/mulan/home/yasheng/Biobank/code/lassosum/lassosum_real2.R
ref=/net/mulan/disk2/yasheng/sample500/frq/merge
val=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
bfilete=/net/mulan/disk2/yasheng/plink_file/genotype/chr
thread=8

for PHENO in 10;do
for cross in 1 3 5 ; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
summ=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/summ/summary_cross${cross}
cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
Rscript ${lassosum} --summ ${summ}.assoc.txt --pheno ${PHENO} --thread ${thread} --cross ${cross}
fi
done
done
