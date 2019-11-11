PHENO=$1
cross=$2
toldpred=/net/mulan/home/yasheng/Biobank/code/ldpred/toldpred.R
py=/net/mulan/home/yasheng/py3/bin/python
ldpred=/net/mulan/home/yasheng/Biobank/program/ldpred/LDpred.py
split=/net/mulan/home/yasheng/Biobank/code/ldpred/split_res2.R
calcr2=/net/mulan/home/yasheng/Biobank/code/ldpred/r2.R
max=/net/mulan/home/yasheng/Biobank/code/ldpred/max.R

let s=200
let ldr=200

## summary 
summ=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/summ/summary_cross${cross}
cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
# Rscript ${toldpred} --gemma ${summ}.assoc.txt  --ldpred ${summ}_LDpred.sumstat

## h2
herit=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/herit/h2_cross${cross}.log
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`
echo herit: ${h2}

## ldpred (cross validation)
n_file=/net/mulan/disk2/yasheng/phenotype_file/v2/pheno_${PHENO}_n_train.txt
n=`sed -n "${cross}, 1p" ${n_file}`
sub=/net/mulan/disk2/yasheng/plink_file/subsample/frq/merge
coord1=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/summary_cv_cross${cross}.HDF5
${py} ${ldpred} coord --gf ${sub} --ssf ${summ}_LDpred.sumstat --out ${coord1} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
ldest=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/ld_cross${cross}_ldr${ldr}
infest=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/esteff_cross${cross}_ldr${ldr}
${py} ${ldpred} gibbs --cf ${coord1} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${h2} --f 0.0001 

## cross validation
pred=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/pheno_cross${cross}
r2=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/r2_cross${cross}
for p in 1.0000e+00 3.0000e-01 1.0000e-01 3.0000e-02 1.0000e-02 3.0000e-03 1.0000e-03 3.0000e-04 1.0000e-04;do
plink-1.9 --bfile ${sub} --score ${infest}_LDpred_p${p}.txt 3 4 7 header sum --out ${pred}_p${p}
Rscript ${calcr2} --pred ${pred}_p${p}.profile --pheno ${PHENO} --r2 ${r2}_p${p}.txt
rm ${pred}_p${p}.log
rm ${pred}_p${p}.nopred
rm ${pred}_p${p}.profile
done
pbestf=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/r2_cross${cross}_pbest.txt
Rscript ${max} --r2 ${r2}_p --pbest ${pbestf} 
pbest=`cat ${pbestf}`
echo pbest: ${pbest}

## delete files
rm ${ldest}*.pkl.gz
rm ${infest}_LDpred*.txt

## ldpred (reference data)
ref2=/net/mulan/disk2/yasheng/sample${s}/frq/merge
coord2=/net/mulan/disk2/yasheng/out_bbj/pheno${PHENO}/ldpred/summary_ref_cross${cross}_ldr${ldr}.HDF5
${py} ${ldpred} coord --gf ${ref2} --ssf ${summ}_LDpred.sumstat --out ${coord2} --N ${n} --ssf-format STANDARD --max-freq-discrep 0.2
${py} ${ldpred} gibbs --cf ${coord2} --ldr ${ldr} --ldf ${ldest} --out ${infest} --N ${n} --h2 ${h2} --f ${pbest}

## split to chr
Rscript ${split} --tot ${infest}_LDpred-inf.txt --sp ${infest}_inf_chr
Rscript ${split} --tot ${infest}_LDpred_p${pbest}.txt --sp ${infest}_pbest_chr

## remove
rm ${coord2}
rm ${ldest}*.pkl.gz
rm ${summ}.assoc.txt
rm ${summ}_LDpred.sumstat

