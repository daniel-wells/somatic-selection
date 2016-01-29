# module load apps/python3/3.3.3/gcc-4.4.7

########################################################################################################
################## Step 1: Calculate total possible nonsynonymous mutations per codon ##################
########################################################################################################

genetic_code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
       "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

# Return effective number of nonsysnoymous sites
def NorS(aa,codon):
    if genetic_code[codon] == aa:
        return 0
    else:
        return 0.33333333333333

# Generate all possible single base substitutions - some redundancy here
from itertools import combinations, product
def generate(s):
    letters = 'ACGT'
    pool = list(s)

    for indices in combinations(range(3), 1):
        for replacements in product(letters, repeat=1):
            skip = False
            for i, a in zip(indices, replacements):
                if pool[i] == a: skip = True
            if skip: continue

            keys = dict(zip(indices, replacements))
            yield ''.join([pool[i] if i not in indices else keys[i] 
                           for i in range(3)])

# Calculate total nonsynonymous sites for a codon
N_per_codon = {}
for codon in genetic_code:
    N_sites = 0
    for new_codon in list(generate(codon)):
        N_sites = N_sites + NorS(aa=genetic_code[codon],codon=new_codon)
    N_per_codon[codon] = N_sites

##################################################################################################################
################## Step 2: Calculate total possible nonsynonymous mutations per coding sequence ##################
##################################################################################################################

#  from Bio.SeqIO.FastaIO https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py
def SimpleFastaParser(handle):
    """Generator function to iterate over Fasta records (as string tuples).
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        # Remove trailing whitespace, and any internal spaces
        # (and any embedded \r which are possible in mangled files
        # when not opened in universal read lines mode)
        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return  # StopIteration

    assert False, "Should not reach this line"


# Generate over records
def fasta_reader(filename):
  with open(filename) as handle:
    for record in SimpleFastaParser(handle):
      yield record

# Calculate total possible nonsynonyous mutations per coding sequence
N_per_cds_total = []
for cds in fasta_reader("data/raw/Homo_sapiens.GRCh37.75.cds.all.fa"):
    N_per_cds = 0
    coding_sequence = cds[1]
    cds_length = len(coding_sequence)
    if cds_length % 3 == 0 and "N" not in coding_sequence:
        for codon in [coding_sequence[i:i+3] for i in range(0, cds_length, 3)]:
            N_per_cds = N_per_cds + N_per_codon[codon]
        N_per_cds = float(N_per_cds)
        error = "None"
    elif cds_length % 3 != 0:
        N_per_cds = "NA"
        error = "Not multiple of 3"
    elif "N" in coding_sequence:
        N_per_cds = "NA"
        error = "N in cds"
    else:
        N_per_cds = "NA"
        error = "Other"
    cds_data = [cds[0],N_per_cds,cds_length,error]
    N_per_cds_total.append(cds_data)

# Write data to tsv file
import csv
with open("data/N_per_transcript.tsv", 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerows(N_per_cds_total)
print "Data written to data/N_per_transcript.tsv"
