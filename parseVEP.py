#!/usr/bin/python
import sys
import vcf
from collections import defaultdict
def fill_annot_dict(annot):
	spl = annot.split("|")
	annot_dict = {}
	annot_dict['transcript'] = spl[2]
	annot_dict['consequence'] = spl[4]
	annot_dict['existing_variation'] = spl[10]
	annot_dict['exon'] = spl[13]
	if spl[25]:
		annot_dict['sift_score'] = float(spl[25][:-1].split("(")[-1])
		annot_dict['sift_cons'] = spl[25].split("(")[0]
	else:
		annot_dict['sift_score'] = "NA"
		annot_dict['sift_cons'] = "NA"
	if spl[26]:
		annot_dict['polyphen_score'] = float(spl[26][:-1].split("(")[-1])
		annot_dict['polyphen_cons'] = spl[26].split("(")[0]
	else:
		annot_dict['polyphen_score'] = "NA"
		annot_dict['polyphen_cons'] = "NA"
	if spl[41]:
		annot_dict['condel_score'] = float(spl[41][:-1].split("(")[-1])
		annot_dict['condel_cons'] = spl[26].split("(")[0]
	else:
		annot_dict['condel_score'] = "NA"
		annot_dict['condel_cons'] = "NA"
	if spl[40]:
		annot_dict['cadd'] = float(spl[40])
	else:
		annot_dict['cadd'] = "NA"
	if spl[45]:
		annot_dict['lof'] = spl[45]
	else:
		annot_dict['lof'] = "NA"
	#var_annot['phyloP100way'] = "NA"
	#var_annot['GERP'] = "NA"
	return annot_dict

	
def get_most_severe_missense(csq):
	high_impact_consequences = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_insertion", "inframe_deletion","missense_variant", "protein_altering_variant"]
	gene_to_annot = {}
	for transcr in csq:
		annot = transcr.split("|")
		
		if annot[4] in high_impact_consequences:
			gene = annot[1] + ";" + annot[23]
			cur_annot_dict = fill_annot_dict(transcr)
			
			if gene in gene_to_annot:
				max_annot = gene_to_annot[gene]
				#print "already exists", prev_annot
				max_condel = max_annot['condel_score']
				if cur_annot_dict['condel_score'] > max_condel:
					gene_to_annot[gene] = cur_annot_dict
			else:
				gene_to_annot[gene] = cur_annot_dict
	
	return gene_to_annot

def summarize_per_gene(var_annot_list):
	sum = defaultdict()
	tools = ['polyphen_cons', 'condel_cons', 'sift_cons', 'lof']
	for t in tools:
		sum[t] = defaultdict(int)
	for var_annot in var_annot_list:
		for t in tools:
			sum[t][var_annot[t]] += 1	
		if (not var_annot['existing_variation']) and (var_annot['condel_cons'] != 'benign'):
			sum['novel_damaging'] += 1

def write_line(gene, var_annot, header, out_file):
	out.write("\t".join(gene.split(";")))
	for col in header:
		out.write("\t" + str(var_annot[col]))
	out.write("\n")
	
def make_table(summary_per_gene, filename):
	header = ['SNP_position', 'existing_variation', 'transcript', 'exon', 'MAF', 'consequence', 'sift_score', 'sift_cons', 'polyphen_score', 'polyphen_cons', 'condel_score', 'condel_cons', 'cadd', 'lof', 'phyloP100way', 'GERP']
	out = open(filename,"w")
	out.write("gene_id\tgene_name\t" + "\t".join(header) + "\n")
	
	for gene, var_annot_list in summary_per_gene.items():
		for var_annot in var_annot_list:
			out.write("\t".join(gene.split(";")))
			for col in header:
				out.write("\t" + str(var_annot[col]))
			out.write("\n")
	out.close()

def add_snp_numbers(maf, snp_nums):
	snp_nums[0] += 1
	if maf <= 0.01:
		snp_nums[1] += 1
	elif maf > 0.05:
		snp_nums[3] += 1
	else:
		snp_nums[2] += 1
	return snp_nums

def add_snps_impact_numbers(var_annot, snp_nums_impact):
	snp_nums_impact[0] += 1
	if var_annot['condel_cons'] == "benign":
		snp_nums_impact[1] += 1
	if var_annot['condel_cons'] == "possibly_damaging":
		snp_nums_impact[2] += 1
	if var_annot['condel_cons'] == "probably_damaging":
		snp_nums_impact[3] += 1

#vcf_reader = vcf.Reader(open("/Volumes/rset1/Botswana/merged_samples/VEP/chr16_test.vcf.gz"))
vcf_reader = vcf.Reader(open(sys.argv[1]))
vep_fields = vcf_reader.infos['CSQ'].desc.split(" ")[-1].split("|")

header = ['SNP_position', 'existing_variation', 'transcript', 'exon', 'MAF', 'consequence', 'sift_score', 'sift_cons', 'polyphen_score', 'polyphen_cons', 'condel_score', 'condel_cons', 'cadd', 'lof', 'phyloP100way', 'GERP']
out_file = open(sys.argv[2],"w")
out.write("gene_id\tgene_name\t" + "\t".join(header) + "\n")

snp_nums = [0,0,0,0] # number of SNPs per MAF bin: overall, <= 0.01, 0.01 < maf <= 0.05, > 0.05
snp_nums_impact = [0, 0, 0] # number of SNPs with high impact, with high impact + damaging; with high impact + damaging + homo alt samples
for record in vcf_reader:
	maf = float(2 * record.num_hom_alt + record.num_het) / (2*(record.num_called + record.num_unknown))
	#cr = float(record.num_called)/ (record.num_called + record.num_unknown)
	snp_nums = add_snp_numbers(maf, snp_nums)
	if 'phyloP100way' in record.INFO:
		cons = record.INFO['phyloP100way'][0]
	else:
		cons = "NA"
	chr_pos = record.CHROM + ":" + str(record.POS)
	gerp = record.INFO['GERP'][0]
	annots = record.INFO['CSQ']
	gene_to_annot = get_most_severe_missense(annots)
	for gene in gene_to_annot:
		var_annot = gene_to_annot[gene]
		var_annot['phyloP100way'] = cons
		var_annot['GERP'] = gerp
		var_annot['MAF'] = str(maf)
		var_annot['SNP_position'] = chr_pos
		#summary_per_gene[gene].append(var_annot)
		write_line(gene, var_annot, header, out_file)
		
out.close()	
#make_table(summary_per_gene, sys.argv[2])

gene, var_annot['SNP_position'], num_high_impact, num_damaging, num_damaging_homo_alt
for gene, var_annot_list in summary_per_gene.items():
	print gene
	for var_annot in var_annot_list:
		print "\t", var_annot
		
rs183360
rs113207608
106588
108554
rs139734214
		
for record in vcf_reader:
	for annot in record	
		if annot.split("|")[42] or annot.split("|")[43] or annot.split("|")[44] or annot.split("|")[45]:
			print "LoF_info: ", annot.split("|")[42], "LoF_flags: ", annot.split("|")[43], "LoF_filter: ", annot.split("|")[44], "LoF", annot.split("|")[45]


0 Allele
1 Gene
2 Feature
3 Feature_type
4 Consequence
5 cDNA_position
6 CDS_position
7 Protein_position
8 Amino_acids
9 Codons
10 Existing_variation
11 AA_MAF
12 EA_MAF
13 EXON
14 INTRON
15 MOTIF_NAME
16 MOTIF_POS
17 HIGH_INF_POS
18 MOTIF_SCORE_CHANGE
19 DISTANCE
20 STRAND
21 CLIN_SIG
22 CANONICAL
23 SYMBOL
24 SYMBOL_SOURCE
25 SIFT
26 PolyPhen
27 GMAF
28 BIOTYPE
29 ENSP
30 DOMAINS
31 CCDS
32 HGVSc
33 HGVSp
34 AFR_MAF
35 AMR_MAF
36 ASN_MAF
37 EUR_MAF
38 PUBMED
39 CADD_RAW
40 CADD_PHRED
41 Condel
42 LoF_info
43 LoF_flags
44 LoF_filter
45 LoF