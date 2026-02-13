###  LDSC
python $ldsc/munge_sumstats.py \
  --sumstats GCST90296696.txt \
  --N N \
  --snp SNP \
  --a1 ALT \
  --a2  REF \
  --p P \
  --frq MAF \
  --signed-sumstats  BETA,0 \
  --merge-alleles w_hm3.snplist \
  --out  GDM_EUR \
  --chunksize 500000 \

python $ldsc/ldsc.py \
 --rg GDM_EUR.sumstats.gz,DAG_EUR.sumstats.gz \
 --ref-ld-chr $ldsc/resource/EUR/baselineLD_v2.3/baselineLD.@ \
 --w-ld-chr $ldsc/resource/EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@ \
 --frqfile-chr $ldsc/resource/EUR/1000G_Phase3_frq/1000G.EUR.QC.@ \
 --out GDM_DAG_EUR
