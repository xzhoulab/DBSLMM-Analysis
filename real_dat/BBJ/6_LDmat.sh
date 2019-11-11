#!/bin/bash

#SBATCH --partition=mulan,nomosix
#SBATCH --time=1-00:00:00
#SBATCH --job-name=LDmat
#SBATCH --mem=2G
#SBATCH --cpus-per-task=3

#SBATCH --array=1-22
#SBATCH --output=/net/mulan/disk2/yasheng/out/LDmat_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/err/LDmat_%a.err
bash
let k=0

LDmat=/net/mulan/home/yasheng/Biobank/code/out_bbj/6_LDmat.R

for((chr=1;chr<23;chr++));do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript ${LDmat} --chr ${chr}
fi
done
