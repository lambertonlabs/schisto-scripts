#need gatk, picardtools, samtools
#schisto variant calling
#assign read groups
java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups I=data/faust/CF3.final.bam O=data/faust/CF3.rg.bam RGLB=CF3 RGPL=illumina RGPU=CF3 RGSM=CF3
java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups I=data/faust/CF7.final.bam O=data/faust/CF7.rg.bam RGLB=CF7 RGPL=illumina RGPU=CF7 RGSM=CF7

#mark duplicates
java -jar picard/build/libs/picard.jar MarkDuplicates I=data/faust/CF7.rg.bam O=data/faust/CF7.rg.md.bam M=data/faust/CF7.M

#index bam files
samtools index /home/ubuntu/data/faust/CF3.rg.md.bam

#index reference genome
samtools faidx /Genomes/Trematoda/schistosoma_mansoni.PRJEA36577.WBPS11.genomic.fa

#call haplotypes
mkdir GVCF

./gatk/gatk HaplotypeCaller \
-R Genomes/Trematoda/schistosoma_mansoni.PRJEA36577.WBPS11.genomic.fa \
-O data/faust/GVCF/CF3.g.vcf \
-I data/faust/CF3.rg.md.bam \
-ERC GVCF

./gatk/gatk HaplotypeCaller \
-R Genomes/Trematoda/schistosoma_mansoni.PRJEA36577.WBPS11.genomic.fa \
-O data/faust/GVCF/CF7.g.vcf \
-I data/faust/CF7.rg.md.bam \
-ERC GVCF
#takes approx 190 minutes per 5X WGS

#consolidate GVCF
mkdir data/faust/GenomicsDB
./gatk/gatk data/faust/GenomicsDBImport  \
---genomicsdb-workspace-path data/faust/GenomicsDB \
-V data/faust/GVCF/CF7.g.vcf \
-V data/faust/GVCF/CF3.g.vcf \
--intervals 20
