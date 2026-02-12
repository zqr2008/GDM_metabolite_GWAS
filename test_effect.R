rm(list = ls())
library(tidyverse)
library(reshape2)
library(ggLD)
filter <- dplyr::filter
select <- dplyr::select


setwd("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/")

mergedoutput <- read.delim("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/join3trimester", header=FALSE, comment.char="#")
colnameofoutput <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ",
  "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P",
  "ERRCODE", "Trait", "Trimester"
)

colnames(mergedoutput) <- colnameofoutput



Supplementary_annofilelipids <- read.csv("C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/349/second_version_manu/KEYSupplementary Table.csv")

annoinfo <- Supplementary_annofilelipids %>%
  dplyr::select(Gene.refGene,SNP,Vname,lipids,lipidname,lipidstype)

filtered <- mergedoutput %>%
  group_by(ID,Trait) %>%
  filter(any(P<5e-8)) %>%
  dplyr::rename(Vname = Trait,SNP=ID) 


Zdiff_sq <- function(betaP, betaM, seP, seM) {
  # delta-method variance of (betaP^2 - betaM^2)
  var_diff <- 4 * (betaP^2 * seP^2 + betaM^2 * seM^2)
  
  Z_D <- (betaP^2 - betaM^2) / sqrt(var_diff)
  P_D <- 2 * (1 - pnorm(abs(Z_D)))
  
  return(list(
    Z = Z_D,
    P = P_D
  ))
}


mergedoutput_spread <- filtered %>%
  mutate(newID= paste0(SNP,"_",Vname)) %>%
  ungroup()%>%
  dplyr::select(newID,BETA,SE,T_STAT,Trimester,P) %>%
  pivot_wider(id_cols = newID,
              names_from =Trimester,
              values_from = c(BETA,SE,T_STAT,P))

mergedoutput_spread_12 <- mergedoutput_spread %>%
  filter(is.na(BETA_1st)==FALSE & is.na(BETA_2nd)==FALSE)

trimester12 <- data.frame()
for (i in (1:nrow(mergedoutput_spread_12))){
  print(mergedoutput_spread[i,1])
  
  zresult <- Zdiff_sq(as.numeric(mergedoutput_spread_12[i,"BETA_1st"]), 
           as.numeric(mergedoutput_spread_12[i,"BETA_2nd"]),
           as.numeric(mergedoutput_spread_12[i,"SE_1st"]),
           as.numeric(mergedoutput_spread_12[i,"SE_2nd"])) 
  
  trimester12[i,1] <- mergedoutput_spread_12[i,1]
  trimester12[i,2] <- zresult$P
  trimester12[i,3] <- "contrast 1st v.s 2nd"
}



mergedoutput_spread_23 <- mergedoutput_spread %>%
  filter(is.na(BETA_2nd)==FALSE & is.na(BETA_3rd)==FALSE)

trimester23 <- data.frame()
for (i in c(1:nrow(mergedoutput_spread_23))){
  print(mergedoutput_spread_23[i,1])
  zresult <- Zdiff_sq(as.numeric(mergedoutput_spread_23[i,"BETA_2nd"]), 
                      as.numeric(mergedoutput_spread_23[i,"BETA_3rd"]),
                      as.numeric(mergedoutput_spread_23[i,"SE_2nd"]),
                      as.numeric(mergedoutput_spread_23[i,"SE_3rd"])) 
  
  trimester23[i,1] <- mergedoutput_spread_23[i,"newID"]
  trimester23[i,2] <- zresult$P
  trimester23[i,3] <- "contrast 2nd v.s 3rd"
}





mergedoutput_spread_13 <- mergedoutput_spread %>%
  filter(is.na(BETA_1st)==FALSE & is.na(BETA_3rd)==FALSE)

trimester13 <- data.frame()
for (i in c(1:nrow(mergedoutput_spread_13))){
  print(mergedoutput_spread_13[i,1])
  zresult <- Zdiff_sq(as.numeric(mergedoutput_spread_13[i,"BETA_1st"]), 
                      as.numeric(mergedoutput_spread_13[i,"BETA_3rd"]),
                      as.numeric(mergedoutput_spread_13[i,"SE_1st"]),
                      as.numeric(mergedoutput_spread_13[i,"SE_3rd"])) 
  
  trimester13[i,1] <- mergedoutput_spread_13[i,"newID"]
  trimester13[i,2] <- zresult$P
  trimester13[i,3] <-"contrast 1st v.s 3rd"
}


trimester12 <- trimester12 %>%
  mutate(FDR=p.adjust(V2,method = "BH"))

trimester23 <- trimester23 %>%
  mutate(FDR=p.adjust(V2,method = "BH"))

trimester13 <- trimester13 %>%
  mutate(FDR=p.adjust(V2,method = "BH"))



allcontrast<- rbind(trimester12,trimester23)
allcontrast<- rbind(allcontrast,trimester13)
