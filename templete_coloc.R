rm(list = ls())
library(coloc)
library(tidyverse)
library(readr)


setwd("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/temporal/")
filelist <- list.files(pattern = "V9.*linear$")


GDMrelated <- read.csv("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/lipdisGMMAT/coloclipid/GDMcoloc/EastAsian_HbA1c_maternal",sep="")



storecolocresult<- data.frame()

x=1
y=0


for (i in filelist){

print(i)
tobecoloc <- read.delim(i,header=FALSE)


colnameofoutput <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ",
  "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P",
  "ERRCODE", "Trait", "Trimester"
)

colnames(tobecoloc) <- colnameofoutput

tobecoloc <-  tobecoloc %>%
mutate(CHROM=as.numeric(CHROM),
       POS= as.numeric(POS),
       SE= as.numeric(SE),
       A1_FREQ = as.numeric(A1_FREQ),
       P=as.numeric(P)) %>%
 dplyr::select(-A1)

subsetsign <- tobecoloc %>%
  dplyr::filter(P<5e-8 & A1_FREQ>0.05)




if(nrow(subsetsign)==0){
 next
}


chrom_bins <- subsetsign %>%
  group_by(CHROM) %>%
  mutate(
    minpos = min(POS),
    maxpos = max(POS),
    sign = if_else(
      (maxpos - minpos) < 5e5,
      "smallestbin",
      paste0(
        "bin",
        floor((POS - minpos) / 5e5) + 1
      )
    )
  ) %>%
  ungroup()


# Select lead SNP per chromosome-bin
chrom_bins <- chrom_bins %>%
  group_by(CHROM, sign) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()



for (j in chrom_bins$ID){

top_B <-  tobecoloc %>%
   filter(ID==j)


# Step 3: Extract region info
target_chr <- top_B$CHROM
target_pos <- top_B$POS  # or BP

# Define window (±500 kb)
window <- 500000

# Step 4: Subset GWAS A to region
gwasA_subset <- GDMrelated %>%
  filter(CHR == target_chr,
         BP >= (target_pos - window),
         BP <= (target_pos + window)) %>%
  dplyr::rename(GDMSE=SE,GDMBETA=BETA)




mergeRSID_NC <- tobecoloc %>%
  mutate(SNP = case_when(str_detect(ID,"chr") ~ paste0(str_split_fixed(ID,":",n=3)[,1],"_",
                                                         str_split_fixed(ID,":",n=3)[,2]),
         TRUE~ ID)) %>%
  inner_join(gwasA_subset,by ="SNP") %>%
  mutate(ALT=tolower(ALT),
         REF=tolower(REF) ) %>%
  filter(A2==ALT & A1==REF )


dataset1 <- list(
  snp = mergeRSID_NC$SNP,
  beta = as.numeric(mergeRSID_NC$GDMBETA),
  varbeta = as.numeric(mergeRSID_NC$GDMSE^2),
  MAF =  as.numeric(mergeRSID_NC$FREQ),
  N = as.numeric(mergeRSID_NC$N),
  type = "quant",
  sdY = 1
)



dataset2 <- list(
  snp =  mergeRSID_NC$SNP,
  beta = as.numeric(mergeRSID_NC$BETA),
  varbeta = as.numeric(mergeRSID_NC$SE^2),
  MAF =  as.numeric(mergeRSID_NC$A1_FREQ),
  N = as.numeric(mergeRSID_NC$OBS_CT),
  type = "quant",
  sdY = 1
)


if(dim(mergeRSID_NC)[1]!=0){
  trimester <- sub(".*_(1st|2nd|3rd)\\..*", "\\1", filelist)
  result <- coloc.abf(dataset1, dataset2)
  tempresult <- as.data.frame(t(as.data.frame(result$summary))) %>%
    mutate(file = i) %>%
    mutate(leadsnp=j) %>%
    mutate(comparetype = "EastAsian_HbA1c_maternal")
  storecolocresult <- rbind(storecolocresult,tempresult)
}


x=x+1


}
}


write.table(storecolocresult,file = paste0("/storeData/mchri/zhaiqiangrong/zhaiqiangrong/GMMAT/lipdisGMMAT/coloclipid/GDMcoloc/result/","V9","EastAsian_HbA1c_maternal_coloceresultforlipid.csv"),
          fileEncoding = "GBK",sep = ",",quote = F)
