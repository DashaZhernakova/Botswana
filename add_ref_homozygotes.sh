#!/bin/bash
set -o pipefail
 f="/Volumes/rset1/Botswana/Analysis/AN75622"
 id=${f##*/}
 vcf=${f}/INDEL-SNP/${id}.SNP.vcf.gz
 echo "Processing $id"
 #mkdir ${f}/INDEL-SNP/per_chr/
 
 for chr in `seq 1 22` "X" "Y" "MT"
 do
  echo "Extracting chr: $chr"
  out_vcf=${f}/INDEL-SNP/per_chr/${id}.SNP.chr${chr}.vcf.gz

  /Users/dzhernakova/Documents/tools/vcftools_0.1.13/bin/vcftools \
        --gzvcf $vcf \
        --chr ${chr} \
        --recode \
        --stdout \
        | bgzip -c > ${out_vcf}
  returnCode=$?
  echo "return code vcftools: $returnCode"      
  
  tabix -p vcf ${out_vcf}

  echo "Getting the complement"
  
  vcf_refhomo=${out_vcf%vcf.gz}ref_homoz.vcf.tmp
  gunzip -c ${out_vcf} > ${vcf_refhomo}
  
  bcftools isec -C \
	/Volumes/rset1/Botswana/merged_samples/all_samples.SNP.PASS.chr${chr}.vcf.gz \
   	${out_vcf} \
   -O v -w 1 -c all \
   | awk 'BEGIN {FS="\t"}; {OFS="\t"}; {if ($1 !~ /^#/) {$6=".";$8="."; $10="0/0:.:.:.:."; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' \
   >> ${vcf_refhomo}
 
  returnCode=$?
  echo "return code bcftools isec: $returnCode"
 
  vcf-sort ${vcf_refhomo} | bgzip -c > ${vcf_refhomo%tmp}
  returnCode=$?
  echo "return code sort: $returnCode"

  if [ $returnCode -eq 0 ]
  then
   echo "removing temp files"
   rm ${vcf_refhomo} 
  fi
done