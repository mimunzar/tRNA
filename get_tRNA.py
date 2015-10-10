#!/usr/bin/env python

# ==========================================================================================
# Desc. : finds tRNA genes inside Escherichia coli CFT073
# In    : genome_sequence.fasta
# Out   : multifasta record of found tRNA genes
# Auth. : Milan Munzar (xmunza00@stud.fit.vutbr.cz)

import sys
import re

# ==========================================================================================
# open fasta file 
# ==========================================================================================

if len(sys.argv) > 1:
  try:
    f = open(sys.argv[1], 'r')
    sequence = f.read()
  except IOError as e:
    print "IOError: {0}".format(e.strerror)
    sys.exit()
  else:
    f.close()
else:
  print "Param error: Missing fasta file!"
  sys.exit(-1)


# ==========================================================================================
# trim the tRNA sequence and create complementary string
# ==========================================================================================

sequence = re.sub(r"(^>.+$)", "", sequence, flags = re.MULTILINE)
sequence = sequence.upper()
sequence = re.sub('\n', '', sequence)
sequence = re.sub('\r', '', sequence)
sequence = re.sub('T', 'U', sequence)

# anti-sense
bp_dict = {'A':'U', 'U' :'A', 'C':'G', 'G':'C'}

c_sequence = ""
for i in range(len(sequence)):
  if bp_dict.__contains__(sequence[i]):
    c_sequence += bp_dict[sequence[i]]
  else:
    c_sequence += sequence[i]
c_sequence = c_sequence[::-1]

# ==========================================================================================
# check for occurences
# ==========================================================================================

pattern = """
[AUCG]{13}                # 1 - 13 
A                         # 14 
(A|G)                     # 15
[AUCG]{1,3}               # 16 - 17A
G                         # 18
[AUCG]{11,14}             # 19 - 31
(A|C|U)                   # 32
U                         # 33
(?P<Anticodon>[AUCG]{3})  # 34 - 36 Anticodon
(A|G)                     # 37
[AUCG]{11,31}             # 38 - 52
GUUC                      # 53 - 56
(G|A)                     # 57
A                         # 58
[AUCG]                    # 59
(U|C)                     # 60
C                         # 61
[AUCG]{12}                # 62 - 73
CCA                       # 74 - 76
"""

re_obj = re.compile(pattern, re.VERBOSE)
tRNAplus = re_obj.finditer(sequence)
tRNAminus = re_obj.finditer(c_sequence)

# ==========================================================================================
# print results
# ==========================================================================================
amino_acids = {
"CGA" : "Ala",
"CGG" : "Ala",
"CGU" : "Ala",
"CGC" : "Ala",
"GCA" : "Arg",
"GCG" : "Arg",
"GCU" : "Arg",
"GCC" : "Arg",
"UCU" : "Arg",
"UCC" : "Arg",
"UUA" : "Asn",
"UUG" : "Asn",
"CUA" : "Asp",
"CUG" : "Asp", 
"ACA" : "Cys",
"ACG" : "Cys",
"CUU" : "Glu",
"CUC" : "Glu",
"GUU" : "Gln",
"GUC" : "Gln",
"CCA" : "Gly",
"CCG" : "Gly",
"CCU" : "Gly",
"CCC" : "Gly",
"GUA" : "His",
"GUG" : "His",
"UAA" : "Ile",
"UAG" : "Ile",
"UAU" : "Ile",
"AAU" : "Leu",
"AAC" : "Leu",
"GAA" : "Leu",
"GAG" : "Leu",
"GAU" : "Leu",
"GAC" : "Leu",
"UUU" : "Lys",
"UUC" : "Lys",
"UAC" : "Met",
"AAA" : "Phe",
"AAG" : "Phe",
"GGA" : "Pro",
"GGG" : "Pro",
"GGU" : "Pro",
"GGC" : "Pro",
"AGA" : "Ser",
"AGG" : "Ser",
"AGU" : "Ser",
"AGC" : "Ser",
"UCA" : "Ser",
"UCG" : "Ser",
"UGA" : "Thr",
"UGG" : "Thr",
"UGU" : "Thr",
"UGC" : "Thr",
"ACC" : "Trp",
"AUA" : "Tyr",
"AUG" : "Tyr",
"CAA" : "Val",
"CAG" : "Val",
"CAU" : "Val",
"CAC" : "Val",
}

# print as multifasta format
i = 1
for match in tRNAplus:
    print ">tRNA_" + str(i) + '|' + amino_acids[match.group("Anticodon")] + match.group("Anticodon") \
    + '|' + str(match.start()) + '|' + str(len(match.group(0))) + '|' + '+' 

    for j in [match.group(0)[k:k + 80] for k in range(0, len(match.group(0)), 80)]:
        print j 

    i = i + 1

for match in tRNAminus:
    print ">tRNA_" + str(i) + '|' + amino_acids[match.group("Anticodon")] + match.group("Anticodon") \
    + '|' + str(len(sequence) - match.end()) + '|' + str(len(match.group(0))) + '|' + '-' 

    for j in [match.group(0)[k:k + 80] for k in range(0, len(match.group(0)), 80)]:
        print j 

    i = i + 1

# EOF
# ==========================================================================================
