#!/usr/bin/env python

# =======================================================================
# Description : compares results from get_tRNA.py with known tRNA genes
# Author      : Milan Munzar (xmunza00@stud.fit.vutbr.cz)
# Params      : FILE found_sequence.multifasta,
#               FILE known_sequence.multifasta
# Output      : list of found headers which are not known, and list of
# known headers which were not found

import sys
import re

# =======================================================================
# load files - found known
# =======================================================================
if len(sys.argv) > 2:
    try:
        f_found = open(sys.argv[1], 'r')
        found = f_found.read()
        f_known = open(sys.argv[2], 'r')
        known = f_known.read()
    except IOError as e:
        print "IOError: {0}".format(e.strerror)
        sys.exit(-1)
    else:
        f_known.close()
        f_found.close()
else:
    print "Param error: Program takes two files as arguments!"
    print "./compare_tRNA.py [found.multifasta] [known.multifasta]"
    sys.exit(-1)

# =======================================================================
# parse headers 
# =======================================================================

found_headers = re.findall("^>(.+)$", found, flags = re.MULTILINE)
parsed_found = []

for i in range(0, len(found_headers)):

    tmp = found_headers[i].split('|')
    p_found = dict().fromkeys(["START", "END", "INDEX"])
    p_found["START"] = int(tmp[2])
    p_found["END"] = int(tmp[2]) + int(tmp[3])
    p_found["INDEX"] = i
    parsed_found.append(p_found)


known_headers = re.findall("^>(.+)$", known, flags = re.MULTILINE)
re_obj = re.compile("(\d+)-(\d+)")
parsed_known = []

for i in range(0, len(known_headers)):

    k_found = dict().fromkeys(["START", "END", "INDEX"])
    match = re_obj.search(known_headers[i])
    k_found["START"] = min(int(match.group(2)), int(match.group(1))) 
    k_found["END"] = max(int(match.group(1)), int(match.group(2))) 
    k_found["INDEX"] = i
    parsed_known.append(k_found)

# =======================================================================
# compare headers, at least 80% overlap
# =======================================================================

forb_found = []
forb_known = []

for i in parsed_found: 
    for j in parsed_known:

        treshold = min(i["END"] - i["START"], j["END"] - j["START"]) * 0.8

        # left overlap
        if i["START"] < j["END"] < i["END"]:
            if j["END"] - i["START"] >= treshold:
                forb_known.append(j["INDEX"]) 
                forb_found.append(i["INDEX"])
            continue

        # right overlap
        if i["START"] < j["START"] < i["END"] :
            if i["END"] - j["START"] >= treshold:
                forb_known.append(j["INDEX"]) 
                forb_found.append(i["INDEX"])
            continue

# =======================================================================
# print output
# =======================================================================

for i in range(len(found_headers)):
    if not forb_found.__contains__(i):
        print found_headers[i]

print ""

for i in range(len(known_headers)):
    if not forb_known.__contains__(i):
        print known_headers[i]

# EOF
# =======================================================================
