template_script <- '
rm(list = ls())
library(GMMAT)
library(tidyverse)


pheno_update <- read.delim("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/pheno_update.txt", header = TRUE)
GRM <- as.matrix(read.table("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/GRM_raw.txt", check.names = FALSE))
pheno_update$ID <- as.character(pheno_update$ID )

pheno <- pheno_update %>% filter(name == "__METAB__")

model2 <- glmmkin(fixed = value ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7+ PC8 + PC9 + PC10 + AGE + BMI + time, data = pheno,
                  kins = GRM, id = "ID", groups = "disease", random.slope = "time",
                  family = gaussian(link = "identity"))

glmm.score(model2, infile = "/storeData/mchri/tangzhuangyuan/tangzhuangyuan/WGS/DPGT/mom/mom_609_snps_hweout_maf5mhc", outfile = "/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/allallall/allgwassum/__METAB___model2_add10pcs_test.txt")
'


# Create scripts
for (metab in levels(factor(pheno_update$name))) {
  this_script <- gsub("__METAB__", metab, template_script)
  file_name <- paste0("run_", metab, ".R")
  writeLines(this_script, file_name)
  cat("Created:", file_name, "\n")
}
