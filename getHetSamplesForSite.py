import vcf
import sys

vcf_reader = vcf.Reader(open(sys.argv[1]))
if sys.argv[2].startswith("rs"):
    rs_id = sys.argv[2]
    for record in vcf_reader:
        if record.ID == rs_id:
            print str(record.ID) + "\t" + str(record.get_hets())
            
else:
    chr=sys.argv[2].split(":")[0]
    pos=sys.argv[2].split(":")[1]

    variants = vcf_reader.fetch(chr, int(pos) - 1, int(pos))

    for var in variants:
        print str(var) + "\t" + str(var.get_hets())
	print [x.sample for x in var.get_hets()]
