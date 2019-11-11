#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=sblupO
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

#SBATCH --array=1-110%10
#SBATCH --output=/net/mulan/disk2/yasheng/out/sblupO_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/sblupO_%a.err

bash
let k=0
let ldr=200
let s=200
let thread=4
gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64
for PHENO in 10 ; do
for cross in 1 2 3 4 5 ; do
for ((chr=1;chr<23;chr++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
herit=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/herit/h2_cross${cross}.log

## heritability
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`

## snp number  
summf=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/summ/summary_cross${cross}
m=`cat ${summf}.assoc.txt | wc -l`
cojo=$(echo "${m}*(1/${h2}-1)" | bc -l)

## est
summ=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/summ/summary_cross${cross}_chr${chr}
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summ}.assoc.txt > ${summ}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summ}.ma
ref=/net/mulan/disk2/yasheng/sample${s}/frq/chr${chr}
est=/net/mulan/disk2/yasheng/out_sample/pheno${PHENO}/sblup/esteff_cross${cross}_chr${chr}
${gcta} --bfile ${ref} --chr ${chr} --cojo-file ${summ}.ma --cojo-sblup ${cojo} --cojo-wind ${ldr} --thread-num ${thread} --out ${est} 

## remove file
rm ${est}.log
rm ${est}.*badsnps
rm ${summ}.ma
fi
done
done
done
