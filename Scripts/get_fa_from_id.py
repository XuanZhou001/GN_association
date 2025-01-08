#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import difflib

"""
根据序列名称获取fasta序列
"""

fa = sys.argv[1] #原始fa文件
result = sys.argv[2] #最终的fa文件
id_file = sys.argv[3] #fa-id文件
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










