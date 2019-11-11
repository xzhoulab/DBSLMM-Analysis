#!/bin/sh

#SBATCH --partition=nomosix,mulan
#SBATCH --time=24:00:00
#SBATCH --job-name=vcf2plink
#SBATCH --mem=5000

#SBATCH --array=1-21
#SBATCH --output=/net/mulan/disk2/yasheng/out/chr%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/chr%a.err
 
bash
let k=0

for ((chr=1;chr<22;chr++)) ; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
VCFFILE=/net/fantasia/home/xzhousph/data/SummaryData/1000GP_Phase3/raw/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
OUTFILE=/net/mulan/disk2/yasheng/EAS/genotype_tmp/chr${chr}
plink-1.9 --vcf ${VCFFILE} --double-id --make-bed --keep-allele-order --out ${OUTFILE}

IDFILE=/net/mulan/disk2/yasheng/EAS/id.txt
SNPFILE=/net/mulan/disk2/yasheng/EAS/snp.txt
INFILE=/net/mulan/disk2/yasheng/EAS/genotype_tmp/chr${chr}
OUTFILE=/net/mulan/disk2/yasheng/EAS/genotype/chr${chr}

plink-1.9 --bfile ${INFILE} --extract ${SNPFILE} --keep ${IDFILE} --maf 0.01 --hwe 0.001 --geno 0 --make-bed --out ${OUTFILE}
fi
done
