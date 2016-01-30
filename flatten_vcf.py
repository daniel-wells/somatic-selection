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