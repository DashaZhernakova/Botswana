#!/usr/bin/python
import sys
import vcf
import copy
import gzip
import argparse
from pyfaidx import Fasta
from itertools import izip
from tqdm import tqdm
    
def swapGenotype(vcf_reader, genotype):
    if (genotype.is_het) or (genotype.gt_type == None): # don't change heterozygous genotypes or missing
        return genotype
    
    new_genotype = copy.copy(genotype)
    gen_data = new_genotype.data
    tmp_dict = gen_data._asdict() # copy genotype contents into a mutable ordered dict
    
    for attr, value in tmp_dict.items():
        if (type(value) is list) and ((len(value) == 2) or (len(value) == 3)): # suppose a list needs to be swapped if has 2 or 3 values
            tmp_value = list(value) #copy to a tmp list
            tmp_value[0] = value[-1]
            tmp_value[-1] = value[0]
            tmp_dict[attr] = tmp_value
    
    if tmp_dict['GT'] == "0/0":
        tmp_dict['GT'] = "1/1"
    elif tmp_dict['GT'] == "1/1":
        tmp_dict['GT'] = "0/0"
    
    new_genotype.data = gen_data._replace(**tmp_dict) # read from dict into genotype call namedtuple
    
    return new_genotype
    
def getRecordWithSwappedGenotypes(vcf_reader, record):
    new_samples = []
    for sample in record.samples:
        new_samples.append(swapGenotype(vcf_reader,sample))
    
    # new record with swapped ref and alt alleles, swapped genotypes and format fields
    return vcf.model._Record(record.CHROM, record.POS, record.ID, record.ALT[0], record.REF, record.QUAL, 
                           record.FILTER, record.INFO, record.FORMAT, 
                           dict([(x,i) for (i,x) in enumerate(vcf_reader.samples)]), 
                           new_samples)
def complement(nucl):
    if nucl == 'A':
        return 'T'
    if nucl == 'T':
        return 'A'
    if nucl == 'G':
        return 'C'
    if nucl == 'C':
        return 'G'
    
def isAT_GCsnp(ref, alt):
    if ref == 'T' and alt == 'A':
        return True
    if ref == 'A' and alt == 'T':
        return True 
    if ref == 'G' and alt == 'C':
        return True
    if ref == 'C' and alt == 'G':
        return True
    return False

def getNumSnps(vcf_name):
    vcf = gzip.open(vcf_name)
    l = vcf.readline()
    num_snps = 1
    while l.startswith("#"):
        l = vcf.readline()
    for l in vcf:
        num_snps += 1
    vcf.close()
    print "Number of SNPs in the input vcf: ", num_snps
    return num_snps
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Swap ref and alt alleles if necessary')
    parser.add_argument('-i','--in', help='VCF input', required=True)
    parser.add_argument('-o','--out', help='VCF output', required=True)
    parser.add_argument('-r','--fasta', help='indexed fasta reference', required=True)
    parser.add_argument('-pb','--progress', action = "store_true", help='show progress bar', default = False)
    args = vars(parser.parse_args())
    in_vcf = args["in"]
    out_vcf = args["out"]
    fasta = args["fasta"]
    num_snps = None
    if "progress" in args:
        num_snps = getNumSnps(in_vcf)
    
    print "input vcf:", in_vcf
    print "output vcf:", out_vcf
    print "reference fasta file:", fasta
    
    #genome = Fasta('/Users/dashazhernakova/Documents/UMCG/hg19/indices/human_g1k_v37.chr10.fa')
    #vcf_reader = vcf.Reader(open("/Users/dashazhernakova/Documents/Doby/Bake-off/test_cmp/chr10.1512KHX-0021_JH-027P_SNP_INDEL.vcf.gz"))
    #vcf_writer = vcf.Writer(open("/Users/dashazhernakova/Documents/Doby/Bake-off/test_cmp/test.vcf.gz", "w"), vcf_reader)
    
    genome = Fasta(fasta)
    vcf_reader = vcf.Reader(open(in_vcf))
    vcf_writer = vcf.Writer(open(out_vcf, "w"),  vcf_reader)
    
    swapped_cnt = 0
    wrong_cnt = 0
    oringinal_cnt = 0
    complement_cnt = 0
    written_cnt = 0
    ct_cnt = 0
    complement_swapped_cnt = 0
    for record in tqdm(vcf_reader, total = num_snps):
        oringinal_cnt += 1
        if record.is_indel: # skip indels
            vcf_writer.write_record(record)
            written_cnt += 1
            continue
        pos = record.start
        chr = record.CHROM
        ref = record.alleles[0]
        if len(record.alleles) > 2: # skip multiallelic sites
            vcf_writer.write_record(record)
            written_cnt += 1
            continue
        alt = record.alleles[1]
        
        if chr not in genome:
            print "WARNING: chromosome", chr, "not in genome, skipping variant"
            continue
        fasta_ref = genome[chr][pos: pos + 1].seq
        
        #correct
        if ref == fasta_ref: # skip correct ref allele cases
            vcf_writer.write_record(record)
            written_cnt += 1
            continue
        
        #skip A/T and G/C SNPs
        if isAT_GCsnp(ref, alt):
            vcf_writer.write_record(record)
            written_cnt += 1
            ct_cnt += 1
            continue
        
        #swap alleles
        if alt == fasta_ref: # need to swap alleles and genotypes
            vcf_writer.write_record(getRecordWithSwappedGenotypes(vcf_reader, record))
            written_cnt += 1
            swapped_cnt += 1
            continue
        
        #complement the alleles
        if ref == complement(fasta_ref):
            record2 = copy.copy(record)
            record2.REF = complement(ref)
            record2.ALT = complement(alt)
            vcf_writer.write_record(record2)
            complement_cnt += 1
            written_cnt += 1
            continue
        
        #complement and swap
        if alt == complement(fasta_ref):
            record2 = copy.copy(record)
            record2.REF = complement(ref)
            record2.ALT = complement(alt)
            
            vcf_writer.write_record(getRecordWithSwappedGenotypes(vcf_reader, record2))
            complement_swapped_cnt += 1
            written_cnt += 1
        else:
            wrong_cnt += 1
            print "WARNING: Something wrong", fasta_ref, record, record.samples
    
    print "Number of variants in the original file: ", oringinal_cnt
    print "Number of variants written to the new file : ", written_cnt
    print "Number of variants swapped ", swapped_cnt
    print "Number of variants needed to be complemented ", complement_cnt
    print "Number of variants needed to be complemented AND swapped ", complement_swapped_cnt
    print "Number of variants with smth wrong (NOT WRITTEN) ", wrong_cnt
    print "Number of C/T SNPs (left intact) ", ct_cnt
    vcf_writer.close()
