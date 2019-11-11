#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simLMM
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-60%10
#SBATCH --output=/net/mulan/disk2/yasheng/out/simLMM_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/simLMM_%a.err

bash
let k=0

let ldr=35
let thread=2

py=/net/mulan/home/yasheng/py3/bin/python
gemma=/net/mulan/home/yasheng/Biobank/program/GEMMA/gemma-0.98-linux-static
gcta=/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64
plink=/usr/cluster/bin/plink-1.9
SBSLMM=/net/mulan/home/yasheng/Biobank/code/slmm2.1/SBSLMM.R
sbslmm=/net/mulan/home/yasheng/Biobank/code/slmm2.1/sbslmm
PLINK=/net/mulan/home/yasheng/plink2_linux_x86_64/plink2
ldsc=/net/mulan/home/yasheng/Biobank/program/ldsc/ldsc.py
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
lassosum=/net/mulan/home/yasheng/Biobank/code/lassosum/lassosum.R

lmmPheno=/net/mulan/home/yasheng/Biobank/code/simulation/lmmPheno3.R
toplinkf=/net/mulan/home/yasheng/Biobank/code/PT/toplinkf.R
clumpf1=/net/mulan/home/yasheng/Biobank/code/PT/toclumpf.R
simpred=/net/mulan/home/yasheng/Biobank/code/PT/simPred.R
max1=/net/mulan/home/yasheng/Biobank/code/PT/max.R
mkldsc=/net/mulan/home/yasheng/Biobank/code/heritability/mkldsc.R
clumpf2=/net/mulan/home/yasheng/Biobank/code/slmm2.1/toclumpf.R
toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
max2=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R

blockf=/net/mulan/home/yasheng/Biobank/data/LDblock/chr1
valIdx=/net/mulan/disk2/yasheng/plink_file/simulation/ref_idx2.txt
trIdx=/net/mulan/disk2/yasheng/plink_file/simulation/train_idx2.txt
teIdx=/net/mulan/disk2/yasheng/plink_file/simulation/test_idx2.txt

for dist in 2 3; do
for herit in 0.1 0.2 0.5 ; do
for ((cross=1;cross<11;cross++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
## simulate genotype data and phenotype data
path=/net/mulan/disk2/yasheng/simulation2/LMM/data/
Rscript ${lmmPheno} --path ${path} --herit ${herit} --cross ${cross} --dist ${dist}

## training and test set
bfile=/net/mulan/disk2/yasheng/simulation2/LMM/data/geno_her${herit}_cross${cross}_dist${dist}
bfiletr=/net/mulan/disk2/yasheng/simulation2/LMM/data/train_her${herit}_cross${cross}_dist${dist}
bfilete=/net/mulan/disk2/yasheng/simulation2/LMM/data/test_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfile} --keep ${trIdx} --make-bed --out ${bfiletr}
plink-1.9 --bfile ${bfile} --keep ${teIdx} --make-bed --out ${bfilete}
rm ${bfiletr}.log
rm ${bfilete}.log

## reference data
bfileval=/net/mulan/disk2/yasheng/simulation2/LMM/data/val_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfile} --keep ${valIdx} --make-bed --out ${bfileval}
rm ${bfileval}.log
rm ${bfile}*

## summary statistics (GEMMA)
cd /net/mulan/disk2/yasheng/simulation2/LMM
summ=summary_her${herit}_cross${cross}_dist${dist}
${gemma} -bfile ${bfiletr} -notsnp -lm 1 -n 1 -o ${summ}
sed -i '1d' /net/mulan/disk2/yasheng/simulation2/LMM/output/${summ}.assoc.txt
rm /net/mulan/disk2/yasheng/simulation2/LMM/output/${summ}.log.txt

## PT 
summf=/net/mulan/disk2/yasheng/simulation2/LMM/output/summary_her${herit}_cross${cross}_dist${dist}
Rscript ${toplinkf} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt
clump1=/net/mulan/disk2/yasheng/simulation2/LMM/PT/summary_her${herit}_cross${cross}_dist${dist}
pred1=/net/mulan/disk2/yasheng/simulation2/LMM/PT/pheno_her${herit}_cross${cross}_dist${dist}
for p in 5e-8 1e-6 1e-4 1e-3 1e-2 5e-2 1e-1 2e-1 5e-1 1.0; do 
Rscript ${clumpf1} --gemma ${summf}.assoc.txt --plink ${summf}.plink.txt --ref ${bfileval} --pth ${p} \
--clump ${clump1}_p${p} --pred ${pred1}
rm ${clump1}_p${p}.clumped
rm ${clump1}_p${p}.log
rm ${pred1}.log
rm ${pred1}.nopred
r21=/net/mulan/disk2/yasheng/simulation2/LMM/PT/r2_her${herit}_cross${cross}_dist${dist}
Rscript ${simpred} --pheno ${pred1}.profile --r2 ${r21}_p${p}.txt
rm ${pred1}.profile
done
pbestf1=/net/mulan/disk2/yasheng/simulation2/LMM/PT/pbest1_her${herit}_cross${cross}_dist${dist}
pbestf2=/net/mulan/disk2/yasheng/simulation2/LMM/PT/pbest2_her${herit}_cross${cross}_dist${dist}
Rscript ${max1} --r2 ${r21}_p --pbest1 ${pbestf1}.txt --pbest2 ${pbestf2}.txt
pbest1=`cat ${pbestf1}.txt`
pbest2=`cat ${pbestf2}.txt`
plink-1.9 --bfile ${bfilete} --score ${clump1}_p${pbest1}.txt 1 2 3 sum --out ${pred1}
mv ${clump1}_p${pbest1}.txt ${clump1}_best.txt
rm ${r21}_p*
rm ${pred1}.log
rm ${pred1}.txt
rm ${pred2}.txt
rm ${clump1}*

m=`cat ${summf}.assoc.txt | wc -l`
n=`cat ${bfiletr}.fam | wc -l`
refld=/net/mulan/disk2/yasheng/sample500/frq/chr1

## SLMM
tmpPath=/net/mulan/disk2/yasheng/simulation2/LMM/slmm/
outPath=/net/mulan/disk2/yasheng/simulation2/LMM/slmm/
Rscript ${SBSLMM} --summary ${summf}.assoc.txt --tmpPath ${tmpPath} --outPath ${outPath} --plink ${plink} --ref ${refld} \
--sbslmm ${sbslmm} --mafMax 0.2 --nsnp ${m} --block ${blockf}.bed --h2 ${herit} --thread ${thread}
est2=/net/mulan/disk2/yasheng/simulation2/LMM/slmm/esteff_summary_her${herit}_cross${cross}_dist${dist}.assoc
pred2=/net/mulan/disk2/yasheng/simulation2/LMM/slmm/pheno_her${herit}_cross${cross}_dist${dist}
plink-1.9 --bfile ${bfilete} --score ${est2}.sbslmm.txt 1 2 4 sum --out ${pred2}
rm ${pred2}.log
rm ${est2}.badsnps

## SBLUP
cojo=$(echo "${m}*(1/${herit}-1)" | bc -l)
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summf}.assoc.txt > ${summf}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summf}.ma
est200=/net/mulan/disk2/yasheng/simulation2/LMM/sblup/esteff_her${herit}_cross${cross}_dist${dist}_200
# est1000=/net/mulan/disk2/yasheng/simulation2/LMM/sblup/esteff_her${herit}_cross${cross}_dist${dist}_1000
${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 200 --thread-num ${thread} --out ${est200}
# ${gcta} --bfile ${refld} --chr 1 --cojo-file ${summf}.ma --cojo-sblup ${cojo} --cojo-wind 1000 --thread-num ${thread} --out ${est1000}
pred3200=/net/mulan/disk2/yasheng/simulation2/LMM/sblup/pheno_her${herit}_cross${cross}_dist${dist}_200
# pred31000=/net/mulan/disk2/yasheng/simulation2/LMM/sblup/pheno_her${herit}_cross${cross}_dist${dist}_1000
plink-1.9 --bfile ${bfilete} --score ${est200}.sblup.cojo 1 2 4 sum --out ${pred3200}
# plink-1.9 --bfile ${bfilete} --score ${est1000}.sblup.cojo 1 2 4 sum --out ${pred31000}
rm ${est200}.log
# rm ${est1000}.log
rm ${est200}*badsnps
# rm ${est1000}*badsnps
# rm ${est200}.sblup.cojo
# rm ${est1000}.sblup.cojo
rm ${pred3200}.log
# rm ${pred31000}.log

## ldpred
Rscript ${toldpred} --gemma ${summf}.assoc.txt  --ldpred ${summf}_LDpred.sumstat
coord1=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/summary_cv_her${herit}_cross${cross}_dist${dist}.HDF5
${py} ${ldpred} coord --gf ${bfileval} --ssf ${summf}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/ld_her${herit}_cross${cross}_dist${dist}
infest=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/esteff_her${herit}_cross${cross}_dist${dist}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} \
--f 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 

pred41=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/pheno_her${herit}_cross${cross}_dist${dist}
r24=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/r2_her${herit}_cross${cross}_dist${dist}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${bfileval} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred41}_p${p}
Rscript ${simpred} --pheno ${pred41}_p${p}.profile --r2 ${r24}_p${p}.txt
rm ${pred41}_p${p}.log
rm ${pred41}_p${p}.nopred
rm ${pred41}_p${p}.profile
done
pbestf4=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/r2_her${herit}_cross${cross}_dist${dist}_best.txt
Rscript ${max2} --r2 ${r24}_p --pbest ${pbestf4} 
pbest4=`cat ${pbestf4}`

coord2=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/summary_ref_her${herit}_cross${cross}_dist${dist}.HDF5
${py} ${ldpred} coord --gf ${refld} --ssf ${summf}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${herit} --f ${pbest4}
predInf4=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/pheno_her${herit}_cross${cross}_dist${dist}_inf
predp4=/net/mulan/disk2/yasheng/simulation2/LMM/ldpred/pheno_her${herit}_cross${cross}_dist${dist}_pbest
plink-1.9 --bfile ${bfilete} --score ${infest}_LDpred-inf.txt 3 4 7 sum --out ${predInf4}
plink-1.9 --bfile ${bfilete} --score ${infest}_LDpred_p${pbest4}.txt 3 4 7 sum --out ${predp4}
rm ${coord1}
rm ${coord2}
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt
rm ${predInf4}.nopred
rm ${predInf4}.log
rm ${predp4}.nopred
rm ${predp4}.log
rm ${r24}*

## lassosum
res5=/net/mulan/disk2/yasheng/simulation2/LMM/lassosum/res_her${herit}_cross${cross}_dist${dist}
Rscript ${lassosum} --summ ${summf}.assoc.txt --ref ${refld} --valid ${bfileval} --test ${bfilete} --res ${res5}

fi
done
done
done
