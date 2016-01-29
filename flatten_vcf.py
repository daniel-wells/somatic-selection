# import req
# r = re.compile(cds["transcript"] + "\|ENSP[0-9]+\|synonymous_variant")
# match = r.search("CONSEQUENCE=AL050302.1|ENSG00000256715|-|AL050302.1-201|ENST00000540061|ENSP00000444145|synonymous_variant|105T>C|D35D;OCCURRENCE=THCA-SA|2|15|0.13333;affected_donors=2;mutation=A>G;project_count=1;tested_donors=9155")
# if match:
#     print match.group(0)
#

# import re
#
#
# for line in fileinput.input():
#     line = re.sub('foo','bar', line.rstrip())
#     [0-9]+\t       367923  MU28611443      C       T       .       .
#     ^.CONSEQUENCE=
#     print(line)
#
#     "ENST[0-9]+\|ENSP[0-9]+\|[a-z]+_variant\|missense_variant|77G>T|G26V


import fileinput
import re

for line in fileinput.input():
    consequences = re.split('CONSEQUENCE=|;OCCURRENCE=', line)
    occourence = consequences[2].strip()
    occourence2 = occourence.replace("|",":")
    occourence3 = occourence2.replace(";","|")
    occourence4 = occourence3.replace("affected_donors=","")
    for transcript in consequences[1].split(','):
        print transcript + '|' + occourence4


#    occourance = occ_to_end.search(line)
#    transcripts2 = transcripts.search(line)
#    print occourance.group(0)
#    print transcripts.group(1)
#    for transcript in transcripts.findall(line):
#        print transcript + "\n"

    #
    # perl -ne 'next unless /ENST/; @l=split /ENST/;
    # foreach $l (@l[1..$#l]){@x=split /\|/, "ENST$l";
    # print join("\t", @x[0..3])."\n"}'
    # simple_somatic_mutation.aggregated.coding.vcf > mutations_by_transcript.tdv