from subprocess import Popen, PIPE

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns."""
    query = '{0}:{1}-{2}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        yield line.strip().split()

import csv
cds_list = csv.DictReader(open("N_per_transcript_clean_head.tsv"),delimiter='\t')
count=0
for cds in cds_list:
    n=0
    s=0
    if cds["chromosome"].startswith("H"):
        print "Error: {0} not a valid chromosome".format(cds["chromosome"])
    else:
        no_variants = True
        for variant in tabix_query("simple_somatic_mutation.aggregated.coding.sort.vcf.gz",cds["chromosome"],cds["start"],cds["end"]):
            no_variants = False
                
            if "missense_variant" in variant[7]:
                n += 1
            elif "synonymous_variant" in variant[7]:
                s += 1
            else:
                print "Error: programme inconsistancy: neither missense or synonymous variant found"
        if no_variants:
            print "Error: No variants found for {0}".format(cds["transcript"])
        else:
            print "{0} {3} Nonsynon:{1} Synon:{2}".format(,n,s,cds["chromosome"])
        count += 1
print count

import req
r = re.compile(cds["transcript"] + "\|ENSP[0-9]+\|synonymous_variant")
match = r.search("CONSEQUENCE=AL050302.1|ENSG00000256715|-|AL050302.1-201|ENST00000540061|ENSP00000444145|synonymous_variant|105T>C|D35D;OCCURRENCE=THCA-SA|2|15|0.13333;affected_donors=2;mutation=A>G;project_count=1;tested_donors=9155")
if match:
    print match.group(0)