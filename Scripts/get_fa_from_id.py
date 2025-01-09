#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import difflib

"""
get fasta from sequence name
"""

fa = sys.argv[1] # input
result = sys.argv[2] # result
id_file = sys.argv[3] # sequence list

W=False
with open(fa, "r") as r, open(result, "w") as w:
    for l in r:
        if l.startswith(">"):
            id = l.strip().lstrip(">").split("\t")[0].split(" ")[0]
            #print(id)
            with open(id_file, "r") as r2:
                for line in r2:
                    if line.strip().split("\t")[0] == id:
                        W = True
                        w.write(l)
                        break
                    else:
                        W = False
        else:
            if W == True:
                w.write(l)
            else:
                pass










