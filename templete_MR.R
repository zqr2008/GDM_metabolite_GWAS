library(TwoSampleMR)
library(data.table)

############################
# USER CONFIGURATION
############################
PLINK2 <- "/storeData/mchri/zhaiqiangrong/zhaiqiangrong/software/plink2"

LD_REF <- "/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/mom_609_snps_hweout"

EXPOSURE_FILE <- "/storeData/mchri/zhaiqiangrong/zhaiqiangrong/temporal/mom_609_snps_hweout_3rd.L535.glm.linear"

OUTCOME_FILE <- "/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/lipdisGMMAT/coloclipid/GDMcoloc/EastAsian_1hGlucose_maternal"

############################
# LIBRARIES
############################
library(TwoSampleMR)
library(data.table)

############################
# 1. READ PLINK2 EXPOSURE GWAS
############################
# PLINK2 glm.linear files contain continuation lines → fread handles this correctly
exp_raw <- as.data.frame(fread(EXPOSURE_FILE))

# Rename #CHROM → CHROM (MANDATORY)
colnames(exp_raw)[colnames(exp_raw) == "#CHROM"] <- "CHROM"

# Keep only additive model
exp_raw <- exp_raw[exp_raw$TEST == "ADD", ]

############################
# 2. FORMAT EXPOSURE FOR TwoSampleMR
############################
exposure_dat <- format_data(
  exp_raw,
  type = "exposure",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "REF",
  eaf_col = "A1_FREQ",
  pval_col = "P",
  samplesize_col = "OBS_CT",
  chr_col = "CHROM",
  pos_col = "POS"
)

exposure_dat$exposure <- "3rd.L535"

# Instrument threshold
exposure_dat <- subset(exposure_dat, pval.exposure < 5e-8)

############################
# 3. REMOVE SNPs NOT IN LD REFERENCE (CRITICAL)
############################
# Read reference BIM
bim <- fread(paste0(LD_REF, ".bim"))
ref_snps <- bim$V2

exposure_dat <- exposure_dat[
  exposure_dat$SNP %in% ref_snps,
]

cat("Exposure SNPs after LD-reference filter:",
    nrow(exposure_dat), "\n")

############################
# 4. LOCAL PLINK2 CLUMPING
############################
clump_input <- exposure_dat[, c("SNP", "pval.exposure")]
colnames(clump_input) <- c("SNP", "P")

write.table(
  clump_input,
  file = "3rd.L535_exposure_for_clumping.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cmd <- paste(
  PLINK2,
  "--bfile", LD_REF,
  "--clump 3rd.L535_exposure_for_clumping.txt",
  "--clump-p1 5e-8",
  "--clump-p2 1",
  "--clump-r2 0.001",
  "--clump-kb 10000",
  "--out 3rd.L535_exposure_clumped"
)

system(cmd)

clumped <- fread(
  "3rd.L535_exposure_clumped.clumps",
  header = TRUE,
  stringsAsFactors = FALSE
)

exposure_dat <- exposure_dat[
  exposure_dat$SNP %in% as.character(clumped$ID),
]


print(clumped)
cat("Exposure SNPs after clumping:",
    nrow(exposure_dat), "\n")

############################
# 5. READ OUTCOME GWAS (GMMAT)
############################
out_raw <- as.data.frame(fread(OUTCOME_FILE))

outcome_dat <- format_data(
  out_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "FREQ",
  pval_col = "PVAL",
  samplesize_col = "N",
  chr_col = "CHR",
  pos_col = "BP"
)

outcome_dat$outcome <- "EastAsian_1hGlucose_maternal"

############################
# 6. HARMONISATION
############################
dat_harmonised <- harmonise_data(
  exposure_dat,
  outcome_dat,
  action = 2
)

cat("SNPs after harmonisation:",
    nrow(dat_harmonised), "\n")

############################
# 7. MENDELIAN RANDOMISATION
############################
mr_results <- mr(dat_harmonised)
print(mr_results)

############################
# 8. SENSITIVITY ANALYSES
############################
print(mr_heterogeneity(dat_harmonised))
print(mr_pleiotropy_test(dat_harmonised))

############################
# 9. INSTRUMENT STRENGTH
############################
dat_harmonised$F_stat <-
  (dat_harmonised$beta.exposure^2) /
  (dat_harmonised$se.exposure^2)

summary(dat_harmonised$F_stat)

############################
# 10. SAVE RESULTS
############################
fwrite(mr_results, "EastAsian_1hGlucose_maternal_3rd.L535.txt", sep = "\t")
