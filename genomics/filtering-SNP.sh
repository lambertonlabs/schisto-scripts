#purge the indels After filtering, kept 997011 out of a possible 1189271 Sites
#filtering freebayes output
vcftools --vcf DP3g95p5maf05.prim.vcf --remove-indels --recode --recode-INFO-all --out SNP.DP3g95p5maf05

# keep variants that have been successfully genotyped in 50% of individuals, min quality  of 30
vcftools --gzvcf raw.vcf.gz --max-missing 0.5 --minQ 30 --recode --recode-INFO-all --out raw.g5mac3
#After filtering, kept 896262 out of a possible 997011 Sites

#depth of 5 (to be retentive), kept 896262 (no diff...?)
vcftools --vcf raw.g5mac3.recode.vcf --minDP 5 --recode --recode-INFO-all --out raw.indg5dp4

#allele balance of between 0.25 to 0.75
~/software/vcflib/bin/vcffilter -f "AB > 0.25 & AB < 0.75 | AB < 0.01" raw.indg5dp4.recode.vcf > indg5dp5ab.vcf &
#kept 825975
