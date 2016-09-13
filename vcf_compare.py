#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 27.03.2016
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
import argparse
from collections import defaultdict

def swapGenotype(gen):
	if gen == "0/0":
		return "1/1"
	if gen == "1/1":
		return "0/0"
	return gen
def swapAlleles(ref1, alt1, ref2, alt2):
	tmp_ref1 = ref1
    tmp_ref2 = ref2
        		
    ref1 = alt1
    ref2 = alt2
    alt1 = tmp_ref1
    alt2 = tmp_ref2
    return ref1, alt1, ref2, alt2

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compare two vcf files.')
    parser.add_argument('-i','--input1', help='VCF input 1', required=True)
    parser.add_argument('-j','--input2', help='VCF input 2', required=True)
    args = vars(parser.parse_args())
    input1 = args["input1"]
    input2 = args["input2"]
    
    with open(input1) as fh:
        data1 = [x.strip().split("\t") for x in fh.readlines() if not x.startswith("#") and x.strip()]
    with open(input2) as fh:
        data2 = [x.strip().split("\t") for x in fh.readlines() if not x.startswith("#") and x.strip()]

    unknown1 = 0
    unknown2 = 0
    rs2data1 = defaultdict(list)
    for item in data1:
        rs = item[2]
        if rs == ".":
            unknown1 += 1
            continue
        ref = item[3]
        alt = item[4]
        meta = item[9].split(":")
        rs2data1[rs] = (ref, alt, meta[0], meta[1:])

    rs2data2 = defaultdict(list)
    for item in data2:
        rs = item[2]
        if rs == ".":
            unknown2 += 1
            continue
        ref = item[3]
        alt = item[4]
        meta = item[9].split(":")
        rs2data2[rs] = (ref, alt, meta[0], meta[1:])

    keys1 = rs2data1.keys()
    keys2 = rs2data2.keys()

    hit = 0
    miss = 0
    gen_same = 0
    gen_miss = 0
    ref_error = 0
    FN = 0
    for key in keys1:
        if key in rs2data2:
            ref1, alt1, gen1, meta1 = rs2data1[key]
            alt1 = alt1.split(",")[0]
            ref2, alt2, gen2, meta2 = rs2data2[key]
            alt2 = alt2.split(",")[0]

            if ref1 != ref2:
            	if (ref1 == alt2) and (ref2 == alt1):
            		swapAlleles(ref1, alt1, ref2, alt2)
            		gen1 = swapGenotype(gen1)
            		gen2 = swapGenotype(gen2)
            		print ref1, alt1, ref2, alt2
            		print gen1, gen2
            		print rs2data1[key]
            		print rs2data2[key]
            		break
            	else:	
                	ref_error += 1
                	continue

            if alt1 == alt2:
                hit += 1
                if gen1 == gen2:
                    gen_same += 1
                else:
                    gen_miss += 1
            else:
                miss += 1
        else:
            FN += 1

    print
    print "D1 size: %s, D2 size %s" % (len(data1), len(data2))
    print "D1 without rs: %s, D2 without rs: %s" % (unknown1, unknown2)
    print
    print "Dataset 1 to Dataset 2"
    print "\tfound correct: %s" % hit
    print "\t\tfound corect genotype: %s" % gen_same
    print "\t\tfound wrong genotype: %s" % gen_miss
    print "\tfound wrong: %s" % miss
    print "\twrong reference: %s" % ref_error
    print
    print "F1", round(2.*hit/(2.*hit+miss+FN),4)
    print "Recall", round(1.*hit/(len(keys1)),4)
    print "Precision", round(1.*gen_same/(hit+miss),4)

    hit = 0
    miss = 0
    gen_same = 0
    gen_miss = 0
    ref_error = 0
    FN = 0
    for key in keys2:
        if key in rs2data1:
            ref1, alt1, gen1, meta1 = rs2data1[key]
            alt1 = alt1.split(",")[0]

            ref2, alt2, gen2, meta2 = rs2data2[key]
            alt2 = alt2.split(",")[0]

            if ref1 != ref2:
            	if (ref1 == alt2) and (ref2 == alt1):
            		tmp_ref1 = ref1
            		tmp_ref2 = ref2
            		
            		ref1 = alt1
            		ref2 = alt2
            		alt1 = tmp_ref1
            		alt2 = tmp_ref2
            	else:
                	ref_error += 1
                	continue

            if alt1 == alt2:
                hit += 1
                if gen1 == gen2:
                    gen_same += 1
                else:
                    gen_miss += 1
            else:
                miss += 1
        else:
            FN += 1

    print
    print "Dataset 1 to Dataset 2"
    print "\tfound correct: %s" % hit
    print "\t\tfound corect genotype: %s" % gen_same
    print "\t\tfound wrong genotype: %s" % gen_miss
    print "\tfound wrong: %s" % miss
    print "\twrong reference: %s" % ref_error
    print
    print "F1", round(2.*hit/(2.*hit+miss+FN),4)
    print "Recall", round(1.*hit/(len(keys2)),4)
    print "Precision", round(1.*gen_same/(hit+miss),4)
