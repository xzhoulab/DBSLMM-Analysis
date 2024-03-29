## SLMM setting
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
  make_option("--p", type="character", default=NULL,
              help="INPUT: proportion of fix effect snp", metavar="numeric"),
  make_option("--prop", type="character", default=NULL,
              help="INPUT: proportion of fix effect heritability", metavar="numeric"),
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
geno_str <- paste0(opt$path, "geno_her", opt$herit, "_cross", opt$cross, 
                   "_dist", opt$dist, "_ps", opt$p, "_prop", opt$prop)
eff_str <- paste0(opt$path, "eff_her", opt$herit, "_cross", opt$cross, 
                  "_dist", opt$dist, "_ps", opt$p, "_prop", opt$prop, ".txt")
pheno_str <- paste0(opt$path, "pheno_her", opt$herit, "_cross", opt$cross,
                    "_dist", opt$dist, "_ps", opt$p, "_prop", opt$prop)


## 
snp_str <- "/net/mulan/disk2/yasheng/plink_file/simulation_tm/snp4.txt"
cp_bed_cmd <- paste0("cp ", snp_info_str, ".bed ", geno_str, ".bed")
cp_bim_cmd <- paste0("cp ", snp_info_str, ".bim ", geno_str, ".bim")
cp_fam_cmd <- paste0("cp ", snp_info_str, ".fam ", geno_str, ".fam")
system(cp_bed_cmd)
system(cp_bim_cmd)
system(cp_fam_cmd)

## fix effect
snp_sel <- data.frame(fread(snp_str))[, 1]
set.seed(seed_str[as.numeric(opt$cross)])
herit_fix <- as.numeric(opt$herit) * as.numeric(opt$prop)
set.seed(seed_str[as.numeric(opt$cross)])
snp_sel_fix_idx <- sample(c(1: length(snp_sel)), floor(length(snp_sel) * as.numeric(opt$p)))
set.seed(seed_str[as.numeric(opt$cross)])
if(opt$dist == "1"){
  eff_fix <- rnorm(length(snp_sel_fix_idx), 0, sqrt(herit_fix / length(snp_sel_fix_idx)))
}
if(opt$dist == "2"){
  eff_fix <- rt(length(snp_sel_fix_idx), 4) / sqrt(length(snp_sel_fix_idx)) * herit_fix
}
if(opt$dist == "3"){
  eff_fix <- rlaplace(length(snp_sel_fix_idx), 0, sqrt(herit_fix/(2*length(snp_sel_fix_idx))))
}
cat ("number of fix snp: ", length(eff_fix), "\n")

## ran effect
herit_ran <- as.numeric(opt$herit) * (1 - as.numeric(opt$prop))
set.seed(seed_str[as.numeric(opt$cross)])
if(opt$dist == "1"){
  eff_ran <- rnorm(length(snp_sel), 0, sqrt(herit_ran / length(snp_sel)))
}
if(opt$dist == "2"){
  eff_ran <- rt(length(snp_sel), 4) / sqrt(length(snp_sel)) * herit_ran
}
if(opt$dist == "3"){
  eff_ran <- rlaplace(length(snp_sel), 0, sqrt(as.numeric(opt$herit)/(2*length(snp_sel))))
}

## 
eff <- vector("numeric", length = length(snp_sel))
eff[snp_sel_fix_idx] <- eff_fix
eff <- eff + eff_ran
write.table(cbind(snp_sel, eff), file = eff_str, col.names = F, row.names = F, quote = F)

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
