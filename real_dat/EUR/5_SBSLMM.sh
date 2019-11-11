#!/bin/bash

#SBATCH --partition=nomosix,mulan
#SBATCH --time=1-00:00:00
#SBATCH --job-name=SBSLMMO
#SBATCH --mem=2G

#SBATCH --array=1-110%10
#SBATCH --output=/net/mulan/disk2/yasheng/out/SBSLMMO_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/SBSLMMO_%a.err
bash
let k=0

plink=/usr/cluster/bin/plink-1.9
SBSLMM=/net/mulan/home/yasheng/Biobank/code/slmm2.1/SBSLMM.R
sbslmm=/net/mulan/home/yasheng/Biobank/code/slmm2.1/sbslmm
thread=4

for PHENO in 10 ; do
for cross in 1 2 3 4 5;do
summf=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/summ/summary_cross${cross}
m=`cat ${summf}.assoc.txt | wc -l`
n_file=/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_${PHENO}_n_train.txt
n=`sed -n "${cross}, 1p" ${n_file}`
herit=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/herit/h2_cross${cross}.log
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`

for ((chr=1;chr<23;chr++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
refld=/net/mulan/disk2/yasheng/sample500/frq/chr${chr}
tmpPath=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/slmm/
outPath=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/slmm/
blockf=/net/mulan/home/yasheng/Biobank/data/LDblock/chr${chr}
Rscript ${SBSLMM} --summary ${summf}_chr${chr}.assoc.txt --tmpPath ${tmpPath} --outPath ${outPath} --plink ${plink} --ref ${refld} \
--sbslmm ${sbslmm} --mafMax 0.2 --nsnp ${m} --block ${blockf}.bed --h2 ${h2} --thread ${thread}
rm ${outPath}esteff_summary_cross${cross}_chr${chr}.*badsnps
fi
done
done
done
