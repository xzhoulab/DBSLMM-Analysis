#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=herit
#SBATCH --mem-per-cpu=2G

#SBATCH --array=1-4
#SBATCH --output=/net/mulan/disk2/yasheng/out/herit_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/herit_%a.err

bash
let k=0

ldsc=/net/mulan/home/yasheng/Biobank/program/ldsc/ldsc.py
mkldsc=/net/mulan/home/yasheng/Biobank/code/heritability/mkldsc.R

source activate ldsc
for PHENO in 1 ; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## summary data for ldsc
summ=/net/mulan/disk2/yasheng/out_sample/raw_data/pheno${PHENO}
ref=/net/mulan/disk2/yasheng/sample500/ldsc/
h2=/net/mulan/disk2/yasheng/out_sample/herit/h2_${PHENO}
## heritability
python ${ldsc} --h2 ${summ}.ldsc.gz --ref-ld-chr ${ref} --w-ld-chr ${ref} --out ${h2}
fi
done
