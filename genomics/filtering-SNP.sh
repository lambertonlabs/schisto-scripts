#FILTER 1
#Mapping quality > 20, only snps not missing data in any individals
#time: 22secs
#variants at start: 1154161
#variants at end: 541015
vcftools --vcf CF-genotyped.vcf --minQ 30 --max-missing 1 --recode --recode-INFO-all --out minQ30max1.vcfll

#FILTER 2
#No INDELß
#time: 16secs
#variants at start: 541015
#variants at end: 491397
vcftools --vcf minQ30max1.vcfll.recode.vcf --remove-indels --recode --recode-INFO-all --out minQ30snp.vcf.recode.vcf

#FILTER 3
#DoC > 5
#time: 16s
#variants at start: 491397
#variants at end: 491397 (?)
vcftools --vcf minQ30snp.vcf.recode.vcf.recode.vcf --minDP 4 --recode --recode-INFO-all --out minQ30dp5snp.vcf



 #haven't finished this yet
#FILTER 4
#filter on allele balances
#time: 16s
#variants at start: 491397
#variants at end: 491397 (?)

./vcflib/bin/vcffilter -s -f "AB > (0.25) & AB < (0.75) | AB < (0.01)" data/faust/variants/unfiltered/minQ30dp5snp.vcf.recode.vcf > minQ30dp5snp.ab.vcf

 bcftools view -g het -i  "AB>0.25 && AB < 0.75" input.vcf.gz

 bcftools view -i 'AB>0.25 && AB < 0.75 | AB < 0.01' minQ30dp5snp.vcf.recode.vcf
