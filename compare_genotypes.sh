#!/bin/bash

#id=TK3016202001505

tail -n+3 $1 | \
while read line
do
id=$line
echo "Processing $id"

chip_vcf=/Volumes/rset1/Botswana_chip/final_with_rel.shared_samples.no_repeats.chr_fixed.fixedAlleles.vcf.gz

time bcftools view -s ${id} \
${chip_vcf} \
-Ou \
| bcftools view \
-e 'GT="0/0" && GT="."' -Oz \
> ${chip_vcf%vcf.gz}nohomoref.${id}.vcf.gz 

tabix -p vcf -f ${chip_vcf%vcf.gz}nohomoref.${id}.vcf.gz 


seq_vcf=/Volumes/rset1/Botswana/snpEff_annotated/botswana_variants_dbsnp.PASS.SNPs.vcf.gz

time bcftools view -s ${id} \
${seq_vcf} \
-Ou \
| bcftools view \
-e 'GT="."' -Oz \
> ${seq_vcf%vcf.gz}${id}.vcf.gz

tabix -f -p vcf ${seq_vcf%vcf.gz}${id}.vcf.gz

time ~/Documents/tools/rtg-tools-3.6.2/rtg vcfeval \
-b ${chip_vcf%vcf.gz}nohomoref.${id}.vcf.gz \
-c ${seq_vcf%vcf.gz}${id}.vcf.gz \
-t ~/Documents/resources/b37/fasta/human_g1k_v37.fa.sdf/ \
-o ./cmp_${id} \
--sample=${id}


echo "Finished, removing vcf files"
rm ${chip_vcf%vcf.gz}nohomoref.${id}.vcf.gz
rm ${seq_vcf%vcf.gz}${id}.vcf.gz

done
