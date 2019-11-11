## LMM setting
library(data.table)
library(plyr)
library(rmutil)
library(optparse)

## parameter setting
args_list = list(
  make_option("--path", type="character", default=NULL,
              help="INPUT: pathway", metavar="numeric"),
  make_option("--herit", type="character", default=NULL,
              help="INPUT: heritability", metavar="numeric"),
  make_option("--cross", type="character", default=NULL,
              help="INPUT: cross", metavar="numeric"), 
  make_option("--dist", type="character", default=NULL,
              help="INPUT: distribution assumption", metavar="numeric")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## 
seed_str <- c(1111, 529, 93, 725, 1988, 2017, 1932, 989, 20170529, 19890725, 
              234, 1908, 23457, 76891, 24357, 876099, 8899977, 65430, 776655, 876445)
snp_info_str <- "/net/mulan/disk2/yasheng/plink_file/simulation/sim_small_sub"

##
geno_str <- paste0(opt$path, "geno_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist)
eff_str <- paste0(opt$path, "eff_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist, ".txt")
pheno_str <- paste0(opt$path, "pheno_her", opt$herit, "_cross", opt$cross, "_dist", opt$dist)

## 
snp_str <- "/net/mulan/disk2/yasheng/plink_file/simulation_tm/snp4.txt"
cp_bed_cmd <- paste0("cp ", snp_info_str, ".bed ", geno_str, ".bed")
cp_bim_cmd <- paste0("cp ", snp_info_str, ".bim ", geno_str, ".bim")
cp_fam_cmd <- paste0("cp ", snp_info_str, ".fam ", geno_str, ".fam")
system(cp_bed_cmd)
system(cp_bim_cmd)
system(cp_fam_cmd)

##
snp_sel <- data.frame(fread(snp_str))[, 1]
set.seed(seed_str[as.numeric(opt$cross)])
if(opt$dist == "1"){
  eff_var <- as.numeric(opt$herit) / length(snp_sel)
  eff <- rnorm(length(snp_sel), 0, sqrt(eff_var))
}
if(opt$dist == "2"){
  eff <- rt(length(snp_sel), 4) / sqrt(length(snp_sel)) * as.numeric(opt$herit)
}
if(opt$dist == "3"){
  eff <- rlaplace(length(snp_sel), 0, sqrt(as.numeric(opt$herit)/(2*length(snp_sel))))
}
eff <- cbind(snp_sel, eff)
write.table(eff, file = eff_str, col.names = F, row.names = F, quote = F)

gcta_str <- "/net/mulan/home/yasheng/summAnnot/software/gcta_1.91.7beta/gcta64"
pheno_cmd <- paste0(gcta_str, " --bfile ", geno_str, " --simu-qt --simu-causal-loci ", eff_str,
                    " --simu-hsq ", opt$herit, " --simu-rep 1 --out ", pheno_str)
system(pheno_cmd)


##
pheno_sim <- data.frame(fread(paste0(pheno_str, ".phen")))[, 3]
fam <- data.frame(fread(paste0(geno_str, ".fam")))
fam[, 6] <- pheno_sim
write.table(fam, paste0(geno_str, ".fam"), col.names = F, row.names = F, quote = F)

##
system(paste0("rm ", pheno_str, ".log"))
system(paste0("rm ", pheno_str, ".par"))
system(paste0("rm ", pheno_str, ".phen"))
system(paste0("rm ", eff_str))
