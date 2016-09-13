#!/usr/bin/python
import sys
import gzip

sampleId = sys.argv[2]
#sampleId = "AV10214"
print "##fileformat=VCFv4.2"
print """##contig=<ID=1,length=2147483645>
##contig=<ID=2,length=2147483645>
##contig=<ID=3,length=2147483645>
##contig=<ID=4,length=2147483645>
##contig=<ID=5,length=2147483645>
##contig=<ID=6,length=2147483645>
##contig=<ID=7,length=2147483645>
##contig=<ID=8,length=2147483645>
##contig=<ID=9,length=2147483645>
##contig=<ID=10,length=2147483645>
##contig=<ID=11,length=2147483645>
##contig=<ID=12,length=2147483645>
##contig=<ID=13,length=2147483645>
##contig=<ID=14,length=2147483645>
##contig=<ID=15,length=2147483645>
##contig=<ID=16,length=2147483645>
##contig=<ID=17,length=2147483645>
##contig=<ID=18,length=2147483645>
##contig=<ID=19,length=2147483645>
##contig=<ID=20,length=2147483645>
##contig=<ID=21,length=2147483645>
##contig=<ID=22,length=2147483645>
##contig=<ID=X,length=2147483645>
##contig=<ID=Y,length=2147483645>
##contig=<ID=MT,length=2147483645>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"""
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sampleId

#with gzip.open("/Users/dashazhernakova/Documents/Doby/Botswana/exome-seq_genotypes/AV10214_SNP_Indel_ANNO.bed.gz") as bed:
with gzip.open(sys.argv[1]) as bed:
	for l in bed:
		spl = l.strip().split("\t")
		ref_AD = str(int(spl[7]) - int(spl[8]))

		genotype = "1/1"
		if spl[5] == "het":
			genotype = "0/1"
				
		print "\t".join([spl[0], str(int(spl[1]) + 1), spl[9], spl[3], spl[4], spl[6], ".", ".", "GT:AD", genotype + ":" + ref_AD + "," + spl[8]])