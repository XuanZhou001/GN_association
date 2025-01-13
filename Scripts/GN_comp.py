#!/usr/bin/env python3

# cite the BiG-SCAPE

from difflib import SequenceMatcher
from scipy.optimize import linear_sum_assignment
from collections import defaultdict
import numpy as np
import os
import sys
import pandas as pd
import subprocess

homo_list = sys.argv[1]
bed_dir = sys.argv[2]
all_pep = sys.argv[3]
all_hmmout = sys.argv[4]
pfam_A_path = sys.argv[5]

def merge_files(file_a, file_b, output_file):
    # merge DSS Index and Jaccard Index into one file
    with open(file_a, 'r') as f:
        lines = f.readlines()
        lines_a = [line.strip().split("\t") for line in lines]
    
    with open(file_b, 'r') as f:
        lines = f.readlines()
        lines_b = [line.strip().split("\t") for line in lines]
    
    if len(lines_a) != len(lines_b):
        raise ValueError("The line numbers are different!")

    num_lines = len(lines_a)
    num_columns = len(lines_a[0])
    
    merged_lines = []
    
    for i in range(num_lines):
        merged_line = []
        for j in range(num_columns):
            if i < j:
                # select DSS index
                merged_line.append(lines_b[i][j])
            else:
                # select Jaccard index
                if j == 0:
                    merged_line.append(lines_a[i][j])
                else:
                    if float(lines_a[i][j]) == 1.0:
                        if float(lines_b[i][j]) != 1.0:
                            merged_line.append(str(0.99999990))
                        else:
                            merged_line.append(lines_a[i][j])
                    else:
                        merged_line.append(lines_a[i][j])

        merged_lines.append(merged_line)

    # the output file
    with open(output_file, 'w') as f:
        for line in merged_lines:
            f.write("\t".join(line) + '\n')


def fasta_parser(handle): # fasta dictionary
    """Parses a fasta file, and stores it in a dictionary.
    Only works if there are no duplicate fasta headers in the fasta file"""
    
    fasta_dict = {}
    header = ""
    with open(handle, 'r') as r:

        for line in r:
            if line[0] == ">":
                header=line.lstrip(">").strip().split()[0]
            else:
                try:
                    fasta_dict[header] += line.strip()
                except KeyError:
                    fasta_dict[header] = line.strip()

    return fasta_dict


def score_expansion(x_string_, y_string_, downstream):
    """
    Input:
    x_string: list of strings. This one tries to expand based on max score
    y_string: reference list. This was already expanded.
    downstream: If true, expansion goes from left to right.
    
    Output:
    max_score, a
    Where a is length of expansion.
    """
    match = 5
    mismatch = -3
    gap = -2
        
    if downstream:
        x_string = x_string_
        y_string = y_string_
    else:
        # Expansion goes upstream. It's easier to flip both slices and proceed
        # as if going downstream.
        x_string = list(reversed(x_string_))
        y_string = list(reversed(y_string_))
    

    # how many gaps to open before calling a mismatch is more convenient?
    max_gaps_before_mismatch = abs(int(match/gap))
    max_score = 0
    score = 0

    pos_y = 0
    a = 0
    b = 0
    for pos_x in range(len(x_string)):
        g = x_string[pos_x]
        
        # try to find g within the rest of the slice
        # This has the obvious problem of what to do if a gene is found _before_
        # the current y_slice (translocation). Duplications could also complicate
        # things
        try:
            match_pos = y_string[pos_y:].index(g)
        except ValueError:
            score += mismatch
        else:
            score += match + match_pos*gap
            pos_y += match_pos + 1 # move pointer one position after match
            
            # Greater or equals. 'equals' because even if max_score isn't 
            # larger, it's an expansion
            if score >= max_score:
                max_score = score
                # keep track of expansion. Account for zero-based numbering
                a = pos_x + 1
    return max_score, a

def cluster_distance_lcs(A, B, A_domlist, B_domlist, dcg_A, dcg_b, core_pos_A, core_pos_b, go_A, go_b, bgc_class, BGCs, AlignedDomainSequences):
    """Compare two clusters using information on their domains, and the 
    sequences of the domains. 
    This version first tries to search for the largest common slices of both BGCs by 
    finding the Longest Common Substring of genes (based on domain content), also
    trying the reverse on one of the BGCs (capital "B" will be used when the 'true'
    orientation of the BGC is found)
    Then the algorithm will try to expand the slices to either side. Once each 
    slice are found, they have to pass one more validation: a core biosynthetic 
    gene (marked in antiSMASH) must be present.
    Capital letters indicate the "true" orientation of the BGC (first BGC, 'A'
    is kept fixed)
    
    Output:
    raw distance, jaccard index, domain sequence similarity, adjacency index, 
    raw DSS non-anchor, raw DSS anchor, number of non-anchor domains in the pair
    (S), number of anchor domains in the pair
    
    """


    mode = 'global'
    bgc_class_weight = {}
    bgc_class_weight["Others"] = (0.01, 0.97, 0.02, 4.0)
    anchor_domains = set()

    Jaccardw, DSSw, AIw, anchorboost = bgc_class_weight["Others"]

    temp_domain_fastas = {}
    
    # Number of genes in each BGC
    lenG_A = len(dcg_A)
    lenG_B = len(dcg_b)
    
    setA = set(A_domlist)
    setB = set(B_domlist)
    intersect = setA & setB
    
    S = 0
    S_anchor = 0
    
    # define the subset of domain sequence tags to include in
    # the DSS calculation. This is done for every domain.
    A_domain_sequence_slice_bottom = defaultdict(int)
    A_domain_sequence_slice_top = defaultdict(int)
    B_domain_sequence_slice_bottom = defaultdict(int)
    B_domain_sequence_slice_top = defaultdict(int)
    
    # Detect totally unrelated pairs from the beginning
    if len(intersect) == 0: # intersection of domains in two GN
        # print("intersect is Zero")
        # Count total number of anchor and non-anchor domain to report in the 
        # network file. Apart from that, these BGCs are totally unrelated.
        for domain in setA: 
            # This is a bit of a hack. If pfam domain ids ever change in size
            # we'd be in trouble. The previous approach was to .split(".")[0]
            # but it's more costly
            # if domain[:7] in anchor_domains: zx
            if domain in anchor_domains:
                S_anchor += len(BGCs[A][domain])
            else:
                S += len(BGCs[A][domain])
                
        for domain in setB:
            # if domain[:7] in anchor_domains: zx
            if domain in anchor_domains:
                S_anchor += len(BGCs[B][domain])
            else:
                S += len(BGCs[B][domain])

        return 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, S, S_anchor, 0, 0, 0, 0


    # initialize domain sequence slices
    # They might change if we manage to find a valid overlap

    for domain in setA:
        A_domain_sequence_slice_bottom[domain] = 0
        A_domain_sequence_slice_top[domain] = len(BGCs[A][domain])
    # print("slice_1")
    # print(A_domain_sequence_slice_top)
    # print(A_domain_sequence_slice_bottom)

    # print("BGCs[B]")
    # print(BGCs[B])

    for domain in setB:
        B_domain_sequence_slice_bottom[domain] = 0
        B_domain_sequence_slice_top[domain] = len(BGCs[B][domain])

    # print(B_domain_sequence_slice_top)
    # print(B_domain_sequence_slice_bottom)


    #initialize domlist borders for AI
    domA_start = 0
    domA_end = len(A_domlist)
    domB_start = 0
    domB_end = len(B_domlist)

    # always find lcs seed to use for offset alignment in visualization
    
    # Compress the list of domains according to gene information. For example:
    # A_domlist = a b c d e f g
    # dcg_a =   1  3  1  2 Number of domains per each gene in the BGC
    # go_a =    1 -1 -1  1 Orientation of each gene
    # A_string = a dcb e fg List of concatenated domains
    # Takes into account gene orientation. This works effectively as putting all
    # genes in the same direction in order to be able to compare their domain content
    A_string = []
    start = 0

    # print("dcg_A")
    # print(dcg_A)
    # print("go_A")
    # print(go_A)
    for g in range(lenG_A):
        domain_count = dcg_A[g]
        if go_A[g] == 1:
            # x[2:] <- small optimization, drop the "PF" from the pfam ids
            # print("".join(x[2:] for x in A_domlist[start:start+domain_count]))
            A_string.append("".join(x[2:] for x in A_domlist[start:start+domain_count]))
        else:
            # print("".join(A_domlist[x][2:] for x in range(start+domain_count-1, start-1 ,-1)))
            A_string.append("".join(A_domlist[x][2:] for x in range(start+domain_count-1, start-1 ,-1)))
        start += domain_count
    # print("A_string:")
    # print(A_string)

    # print("dcg_b")
    # print(dcg_b)
    # print("go_b")
    # print(go_b)
    b_string = []
    start = 0
    for g in range(lenG_B):
        domain_count = dcg_b[g]
        if go_b[g] == 1:
            b_string.append("".join(x[2:] for x in B_domlist[start:start+domain_count]))
        else:
            b_string.append("".join(B_domlist[x][2:] for x in range(start+domain_count-1, start-1 ,-1)))
        start += domain_count
    # print("b_string:")
    # print(b_string)
        
    b_string_reverse = list(reversed(b_string))
        
    seqmatch = SequenceMatcher(None, A_string, b_string)
    # a: start position in A
    # b: start position in B
    # s: length of the match
    a, b, s = seqmatch.find_longest_match(0, lenG_A, 0, lenG_B)
    #print(A, B)
    # print(a, b, s)
    
    seqmatch = SequenceMatcher(None, A_string, b_string_reverse)
    ar, br, sr = seqmatch.find_longest_match(0, lenG_A, 0, lenG_B)
    # print(ar, br, sr)
    
    # We need to keep working with the correct orientation
    if s > sr or (s == sr and go_A[a] == go_b[b]):
        # print("1")
        dcg_B = dcg_b
        B_string = b_string
        # note: these slices are in terms of genes, not domains (which are 
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        sliceStartA = a
        sliceStartB = b
        sliceLengthA = s
        sliceLengthB = s
        
        reverse = False
        b_name = B
        
    elif s < sr:
        # print("2")
        dcg_B = list(reversed(dcg_b))
        B_string = b_string_reverse
        
        sliceStartA = ar
        sliceStartB = br
        sliceLengthA = sr
        sliceLengthB = sr
        
        # We'll need to know if we're working in reverse so that the start 
        # postion of the final slice can be transformed to the original orientation
        reverse = True
        b_name = B + "*"
        
    # special cases: s == sr
    
    # if only one gene matches, choose the one with the most domains
    elif s == 1:
        # print("3")
        seqmatch = SequenceMatcher(None, A_string, b_string)
        max_domains = 0
        x = 0   # index in A with the gene with most domains
        y = 0   # index in B with the gene with most domains
        for a, b, z in seqmatch.get_matching_blocks():
            if z != 0:
                # if DomainCountGene[A][a] > max_domains:
                if dcg_A[a] > max_domains: # zx change
                    x = a
                    y = b
                    # max_domains = DomainCountGene[A][a]
                    max_domains = dcg_A[a] # zx change
        
        # note: these slices are in terms of genes, not domains (which are 
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        sliceStartA = x
        sliceLengthA = 1
        sliceLengthB = 1
        
        # if BGCGeneOrientation[A][x] == BGCGeneOrientation[B][y]: # zx
        if go_A[x] == go_b[y]: #zx    
            sliceStartB = y
            dcg_B = dcg_b
            B_string = b_string
            reverse = False
            b_name = B
        else:
            # sliceStartB = len(BGCGeneOrientation[B]) - y - 1
            sliceStartB = len(go_b) - y - 1 # zx change
            dcg_B = list(reversed(dcg_b))
            B_string = b_string_reverse
            reverse = True
            b_name = B + "*"
               
    # s == sr and (s == 0 or s > 1)
    else:
        # print("4")
        dcg_B = dcg_b
        B_string = b_string
        # note: these slices are in terms of genes, not domains (which are 
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        sliceStartA = a
        sliceStartB = b
        sliceLengthA = s
        sliceLengthB = s
        
        reverse = False
        b_name = B
        
        
    lcsStartA = sliceStartA
    lcsStartB = sliceStartB
    seedLength = sliceLengthA
    # 以上三个参数为最终结果的一部分
    
    
    # paint stuff on screen
    #if sliceStartB > sliceStartA:
        #offset_A = sliceStartB - sliceStartA
        #offset_B = 0
    #else:
        #offset_A = 0
        #offset_B = sliceStartA - sliceStartB
    #print(A, B)
    #print("  "*offset_A + " ".join(map(str,dcg_A[:sliceStartA])) + "[" + " ".join(map(str,dcg_A[sliceStartA:sliceStartA+sliceLengthA])) + "]" + " ".join(map(str,dcg_A[sliceStartA+sliceLengthA:])) + "\t" + A )
    #print("  "*offset_B + " ".join(map(str,dcg_B[:sliceStartB])) + "[" + " ".join(map(str,dcg_B[sliceStartB:sliceStartB+sliceLengthB])) + "]" + " ".join(map(str,dcg_B[sliceStartB+sliceLengthB:])) + "\t" + b_name)
    ##print(sliceStartA, sliceStartB, sliceLengthA)
    #print("")

    if mode=="glocal" or (mode=="auto" and (bgc_info[A].contig_edge or bgc_info[B].contig_edge)):
        # print("mode")
        #X: bgc that drive expansion
        #Y: the other bgc
        # forward: True if expansion is to the right
        # returns max_score, final positions for X and Y
        #score_expansion(X_string, dcg_X, Y_string, dcg_Y, downstream=True/False)
            
        # Expansion is relatively costly. We ask for a minimum of 3 genes
        # for the core overlap before proceeding with expansion.
        biosynthetic_hit_A = False
        for biosynthetic_position in core_pos_A:
            # print(biosynthetic_position) #zx
            if biosynthetic_position >= sliceStartA and biosynthetic_position <= (sliceStartA+sliceLengthA):
                biosynthetic_hit_A = True
                break
        if sliceLengthA >= 3 or biosynthetic_hit_A:
            # LEFT SIDE
            # Find which bgc has the least genes to the left. If both have the same 
            # number, find the one that drives the expansion with highest possible score
            if sliceStartA == sliceStartB:
                # assume complete expansion of A, try to expand B
                score_B, sbB = score_expansion(B_string[:sliceStartB], A_string[:sliceStartA], False)
                # assume complete expansion of B, try to expand A
                score_A, saA = score_expansion(A_string[:sliceStartA], B_string[:sliceStartB], False)
                
                if score_A > score_B or (score_A == score_B and saA > sbB):
                    sliceLengthA += saA
                    sliceLengthB += len(B_string[:sliceStartB])
                    sliceStartA -= saA
                    sliceStartB = 0
                else:
                    sliceLengthA += len(A_string[:sliceStartA])
                    sliceLengthB += sbB
                    sliceStartA = 0
                    sliceStartB -= sbB

            else:
                # A is shorter upstream. Assume complete extension. Find B's extension
                if sliceStartA < sliceStartB:
                    score_B, sb = score_expansion(B_string[:sliceStartB], A_string[:sliceStartA], False)
                    
                    sliceLengthA += len(A_string[:sliceStartA])
                    sliceLengthB += sb
                    sliceStartA = 0
                    sliceStartB -= sb
                else:
                    score_A, sa = score_expansion(A_string[:sliceStartA], B_string[:sliceStartB], False)
                    
                    sliceLengthA += sa
                    sliceLengthB += len(B_string[:sliceStartB])
                    sliceStartA -= sa
                    sliceStartB = 0

            # RIGHT SIDE
            # check which side is the shortest downstream. If both BGCs have the same
            # length left, choose the one with the best expansion
            downstream_A = lenG_A - sliceStartA - sliceLengthA
            downstream_B = lenG_B - sliceStartB - sliceLengthB
            if downstream_A == downstream_B:
                # assume complete extension of A, try to expand B
                score_B, xb = score_expansion(B_string[sliceStartB+sliceLengthB:], A_string[sliceStartA+sliceLengthA:], True)
                # assume complete extension of B, try to expand A
                score_A, xa = score_expansion(A_string[sliceStartA+sliceLengthA:], B_string[sliceStartB+sliceLengthB:], True)
                
            
                if (score_A == score_B and xa > xb) or score_A > score_B:
                    sliceLengthA += xa
                    sliceLengthB += len(B_string[sliceStartB+sliceLengthB:])
                else:
                    sliceLengthA += len(A_string[sliceStartA+sliceLengthA:])
                    sliceLengthB += xb
                    
            else:
                if downstream_A < downstream_B:
                    # extend all of remaining A
                    score_B, xb = score_expansion(B_string[sliceStartB+sliceLengthB:], A_string[sliceStartA+sliceLengthA:], True)
                    
                    sliceLengthA += len(A_string[sliceStartA+sliceLengthA:])
                    sliceLengthB += xb
                    
                else:
                    score_A, xa = score_expansion(A_string[sliceStartA+sliceLengthA:], B_string[sliceStartB+sliceLengthB:], True)
                    sliceLengthA += xa
                    sliceLengthB += len(B_string[sliceStartB+sliceLengthB:])
        
        #print("  "*offset_A + " ".join(map(str,dcg_A[:sliceStartA])) + "[" + " ".join(map(str,dcg_A[sliceStartA:sliceStartA+sliceLengthA])) + "]" + " ".join(map(str,dcg_A[sliceStartA+sliceLengthA:])) + "\t" + A )
        #print("  "*offset_B + " ".join(map(str,dcg_B[:sliceStartB])) + "[" + " ".join(map(str,dcg_B[sliceStartB:sliceStartB+sliceLengthB])) + "]" + " ".join(map(str,dcg_B[sliceStartB+sliceLengthB:])) + "\t" + b_name)

        #print()
        #print(A_string)
        #print(B_string)
        #print()
        if min(sliceLengthA, sliceLengthB) >= 5:
            # First test passed. Find if there is a biosynthetic gene in both slices
            # (note that even if they are, currently we don't check whether it's 
            # actually the _same_ gene)
            biosynthetic_hit_A = False
            biosynthetic_hit_B = False
            
            for biosynthetic_position in core_pos_A:
                if biosynthetic_position >= sliceStartA and biosynthetic_position <= (sliceStartA+sliceLengthA):
                    biosynthetic_hit_A = True
                    break
            
            # return to original orientation if needed
            if reverse:
                sliceStartB = lenG_B - sliceStartB - sliceLengthB
                
            # using original core_pos_b
            for biosynthetic_position in core_pos_b:
                if biosynthetic_position >= sliceStartB and biosynthetic_position <= (sliceStartB + sliceLengthB):
                    biosynthetic_hit_B = True
                    break
                
            # finally...
            if biosynthetic_hit_A and biosynthetic_hit_B:
                domA_start = sum(dcg_A[:sliceStartA])
                domA_end = domA_start + sum(dcg_A[sliceStartA:sliceStartA+sliceLengthA])
                setA = set(A_domlist[domA_start:domA_end])
                
                domB_start = sum(dcg_b[:sliceStartB])
                domB_end = domB_start + sum(dcg_b[sliceStartB:sliceStartB+sliceLengthB])
                setB = set(B_domlist[domB_start:domB_end])
                
                intersect = setA & setB
                
                # re-adjust the indices for each domain so we get only the sequence
                # tags in the selected slice. First step: find out which is the 
                # first copy of each domain we're using
                for domain in A_domlist[:domA_start]:
                    A_domain_sequence_slice_bottom[domain] += 1
                    # print(domain)
                    # print(A_domain_sequence_slice_bottom[domain])
                for domain in B_domlist[:domB_start]:
                    B_domain_sequence_slice_bottom[domain] += 1
                    # print(domain)
                    # print(B_domain_sequence_slice_bottom[domain])
                    
                # Step 2: work with the last copy of each domain. 
                # Step 2a: make top = bottom
                for domain in setA:
                    A_domain_sequence_slice_top[domain] = A_domain_sequence_slice_bottom[domain]
                    # print(domain)
                    # print(A_domain_sequence_slice_top[domain])
                for domain in setB:
                    B_domain_sequence_slice_top[domain] = B_domain_sequence_slice_bottom[domain]
                    # print(domain)
                    # print(B_domain_sequence_slice_top[domain])
                
                # Step 2b: increase top with the domains in the slice
                for domain in A_domlist[domA_start:domA_end]:
                    A_domain_sequence_slice_top[domain] += 1
                    # print(domain)
                    # print(A_domain_sequence_slice_top[domain])
                for domain in B_domlist[domB_start:domB_end]:
                    B_domain_sequence_slice_top[domain] += 1
                    # print(domain)
                    # print(B_domain_sequence_slice_top[domain])
                    
            #else:
                #print(" - - Not a valid overlap found - - (no biosynthetic genes)\n")
                
        #else:
                #print(" - - Not a valid overlap found - - (shortest slice not large enough)\n")
                    
    # JACCARD INDEX
    print(intersect)
    print(setA | setB)
    Jaccard = len(intersect) / len(setA | setB)
    print(Jaccard)


    # DSS INDEX
    #domain_difference: Difference in sequence per domain. If one cluster does
    # not have a domain at all, but the other does, this is a (complete) 
    # difference in sequence 1. If both clusters contain the domain once, and 
    # the sequence is the same, there is a seq diff of 0.
    #S: Max occurence of each domain
    domain_difference_anchor,S_anchor = 0,0
    domain_difference,S = 0,0
        
    not_intersect = setA.symmetric_difference(setB)
    # print("slice_2")
    # print(A_domain_sequence_slice_top)
    # print(A_domain_sequence_slice_bottom)
    # print(B_domain_sequence_slice_top)
    # print(B_domain_sequence_slice_bottom)

        
    # Case 1
    #no need to look at seq identity, since these domains are unshared
    # print("num_unshared")
    for unshared_domain in not_intersect:
        # print(unshared_domain)
        # print(A_domain_sequence_slice_top)
        # print(B_domain_sequence_slice_top)
        # print(A_domain_sequence_slice_bottom)
        # print(B_domain_sequence_slice_bottom)
        #for each occurence of an unshared domain do domain_difference += count 
        # of domain and S += count of domain
        unshared_occurrences = []
        # try:
        #     print(A_domain_sequence_slice_top[unshared_domain])
        #     print(A_domain_sequence_slice_bottom[unshared_domain])
        # except KeyError:
        #     print(B_domain_sequence_slice_top[unshared_domain])
        #     print(B_domain_sequence_slice_bottom[unshared_domain])

        # try:
        #     num_unshared = A_domain_sequence_slice_top[unshared_domain] - A_domain_sequence_slice_bottom[unshared_domain]
        # except KeyError:
        #     num_unshared = B_domain_sequence_slice_top[unshared_domain] - B_domain_sequence_slice_bottom[unshared_domain]
        if unshared_domain in A_domain_sequence_slice_top:
            num_unshared = A_domain_sequence_slice_top[unshared_domain] - A_domain_sequence_slice_bottom[unshared_domain]
        else:
            num_unshared = B_domain_sequence_slice_top[unshared_domain] - B_domain_sequence_slice_bottom[unshared_domain]
        # print(num_unshared)    
        # don't look at domain version, hence the split
        # if unshared_domain[:7] in anchor_domains: zx
        if unshared_domain in anchor_domains: #zx
            domain_difference_anchor += num_unshared
        else:
            domain_difference += num_unshared
                    
    S = domain_difference # can be done because it's the first use of these
    S_anchor = domain_difference_anchor
    # print("domain_difference:")
    # print(S)
    # print("domain_difference_anchor")
    # print(S_anchor)
        
    # Cases 2 and 3 (now merged)
    for shared_domain in intersect:
        # print("shared_domain")
        # print(shared_domain)
        specific_domain_list_A = BGCs[A][shared_domain]
        specific_domain_list_B = BGCs[B][shared_domain]
        
        num_copies_a = A_domain_sequence_slice_top[shared_domain] - A_domain_sequence_slice_bottom[shared_domain]
        num_copies_b = B_domain_sequence_slice_top[shared_domain] - B_domain_sequence_slice_bottom[shared_domain]
        
        temp_domain_fastas.clear()
        
        accumulated_distance = 0
            
        # Fill distance matrix between domain's A and B versions
        DistanceMatrix = np.ndarray((num_copies_a,num_copies_b))
        # print(specific_domain_list_A) # ZX
        # print(specific_domain_list_B) # ZX
        for domsa in range(num_copies_a):
            for domsb in range(num_copies_b):
                sequence_tag_a = specific_domain_list_A[domsa + A_domain_sequence_slice_bottom[shared_domain]]
                sequence_tag_b = specific_domain_list_B[domsb + B_domain_sequence_slice_bottom[shared_domain]]
                
                seq_length = 0
                matches = 0
                gaps = 0
                
                # print(sequence_tag_a + "  " + sequence_tag_b) #ZX
                aligned_seqA = AlignedDomainSequences[sequence_tag_a]
                aligned_seqB = AlignedDomainSequences[sequence_tag_b]
                    
                # - Calculate aligned domain sequences similarity -
                # Sequences *should* be of the same length unless something went
                # wrong elsewhere
                if len(aligned_seqA) != len(aligned_seqB):
                    print("\tWARNING: mismatch in sequences' lengths while calculating sequence identity ({})".format(shared_domain))
                    print("\t  Specific domain 1: {} len: {}".format(sequence_tag_a, str(len(aligned_seqA))))
                    print("\t  Specific domain 2: {} len: {}".format(sequence_tag_b, str(len(aligned_seqB))))
                    seq_length = min(len(aligned_seqA), len(aligned_seqB))
                else:
                    seq_length = len(aligned_seqA)
                    
                for position in range(seq_length):
                    if aligned_seqA[position] == aligned_seqB[position]:
                        if aligned_seqA[position] != "-":
                            matches += 1
                        else:
                            gaps += 1
                            
                DistanceMatrix[domsa][domsb] = 1 - ( matches/(seq_length-gaps) )
                
        #Only use the best scoring pairs
        BestIndexes = linear_sum_assignment(DistanceMatrix)
        accumulated_distance = DistanceMatrix[BestIndexes].sum()
        
        # the difference in number of domains accounts for the "lost" (or not duplicated) domains
        sum_seq_dist = (abs(num_copies_a-num_copies_b) + accumulated_distance)  #essentially 1-sim
        normalization_element = max(num_copies_a, num_copies_b)

           
        # if shared_domain[:7] in anchor_domains: #zx
        if shared_domain in anchor_domains: #zx
            S_anchor += normalization_element
            domain_difference_anchor += sum_seq_dist
        else:
            S += normalization_element
            domain_difference += sum_seq_dist
    # print("S_anchor, S:")
    # print(S_anchor, S)
    # print("domain_difference, domain_difference_anchor:")
    # print(domain_difference, domain_difference_anchor)

    if S_anchor != 0 and S != 0: # 决定DSS的最后一个循环
        DSS_non_anchor = domain_difference / S
        DSS_anchor = domain_difference_anchor / S_anchor
        
        # Calculate proper, proportional weight to each kind of domain
        non_anchor_prct = S / (S + S_anchor)
        anchor_prct = S_anchor / (S + S_anchor)
        
        # boost anchor subcomponent and re-normalize
        non_anchor_weight = non_anchor_prct / (anchor_prct*anchorboost + non_anchor_prct)
        anchor_weight = anchor_prct*anchorboost / (anchor_prct*anchorboost + non_anchor_prct)

        # Use anchorboost parameter to boost percieved rDSS_anchor
        DSS = (non_anchor_weight*DSS_non_anchor) + (anchor_weight*DSS_anchor)
        
    elif S_anchor == 0:
        DSS_non_anchor = domain_difference / S
        DSS_anchor = 0.0
        
        DSS = DSS_non_anchor
        
    else: #only anchor domains were found
        DSS_non_anchor = 0.0
        DSS_anchor = domain_difference_anchor / S_anchor
        
        DSS = DSS_anchor
 
    DSS = 1-DSS #transform into similarity
    # print("DSS:")
    # print(DSS)
 

    # ADJACENCY INDEX
    # calculates the Tanimoto similarity of pairs of adjacent domains
    
    if len(A_domlist[domA_start:domA_end]) < 2 or len(B_domlist[domB_start:domB_end]) < 2:
        AI = 0.0
    else:
        setA_pairs = set()
        for l in range(domA_start, domA_end-1):
            setA_pairs.add(tuple(sorted([A_domlist[l],A_domlist[l+1]])))
        
        setB_pairs = set()
        for l in range(domB_start, domB_end-1):
            setB_pairs.add(tuple(sorted([B_domlist[l],B_domlist[l+1]])))

        # same treatment as in Jaccard
        AI = len(setA_pairs & setB_pairs) / len(setA_pairs | setB_pairs)

    Distance = 1 - (Jaccardw * Jaccard) - (DSSw * DSS) - (AIw * AI)
    
    # This could happen due to numerical innacuracies
    if Distance < 0.0:
        if Distance < -0.000001: # this definitely is something else...
            print("Negative distance detected!")
            print(Distance)
            print("{} - {}".format(A, B))
            print("J: {}\tDSS: {}\tAI: {}".format(str(Jaccard), str(DSS), str(AI)))
            print("Jw: {}\tDSSw: {}\tAIw: {}".format(str(Jaccardw), str(DSSw), str(AIw)))
        Distance = 0.0
        
    rev = 0.0
    if reverse:
        rev = 1.0
        
    return Distance, Jaccard, DSS, AI, DSS_non_anchor, DSS_anchor, S, S_anchor, lcsStartA, lcsStartB, seedLength, rev



def get_fa_from_id(fa, result, g_list): # 修改为输入list
    # 根据id提取蛋白序列
    W=False
    with open(fa, "r") as r, open(result, "w") as w:
        for l in r:
            if l.startswith(">"):
                Id = l.strip().lstrip(">").split("\t")[0].split(" ")[0]
                #print(id)
                if Id in g_list:
                    W = True
                    w.write(l)
                    # break
                else:
                    W = False
                # with open(id_file, "r") as r2:
                #     for line in r2:
                #         if line.strip().split("\t")[0] == id:
                #             W = True
                #             w.write(l)
                #             break
                #         else:
                #             W = False
            else:
                if W == True:
                    w.write(l)
                else:
                    pass

def stockholm_parser(stkFile):
    reference = ""
    algnDict = defaultdict(str)
    algnFile = stkFile[:-3] + 'algn'
    if not os.path.isfile(algnFile): # prevents overwriting algn files
        with open(stkFile, 'r') as infile:
            for l in infile:
                line = l.strip()
                if line.startswith("#=GC RF"):
                    reference += line[7:].strip()
                elif line == "":
                    continue
                elif line[0] == "/" or line[0] == "#":
                    continue
                else:
                    a = line.split(" ")
                    header = a[0]
                    algn = a[-1]
                    algnDict[header] += algn
                    
        # get start-end coordinates of every "x" island (original consensus)
        # in the reference
        state_reference = False
        slicing_tuples = []
        for pos in range(len(reference)):
            if reference[pos] == "x" and not state_reference:
                state_reference = True
                start = pos
            if reference[pos] == "." and state_reference:
                state_reference = False
                slicing_tuples.append((start,pos))
        if state_reference:
            slicing_tuples.append((start, len(reference)))
    
    if len(algnDict) > 0:
        with open(algnFile, "w") as outfile:
            for header in algnDict:
                sequence = ""
                for a,b in slicing_tuples:
                    sequence += algnDict[header][a:b]
                outfile.write(">{}\n".format(header))
                outfile.write(sequence + "\n")
    return


def DSS_and_JaccardIndex(homo_list, bed_dir, all_pep, all_hmmout, pfam_A_path):
    domout_DB = "domout"
    domains_seq_folder = "domain_seq"
    if os.path.exists(domout_DB):
        pass
    else:
        os.mkdir(domout_DB)

    if os.path.exists(domains_seq_folder):
        pass
    else:
        os.mkdir(domains_seq_folder)

    # 存储所有pfam id
    DomainList = {} 
    # 记录每个基因有几个pfam结果
    DomainNumberGene = {}
    # 记录核心基因的index
    CoreGeneIndex = {}
    # 记录基因的方向
    GeneOrientation = {}
    # 以PFAMid为键的值
    BGCs = {}


    Homo_list = []
    domain_list = []

    gene_sp_dic = {}


    with open(homo_list, 'r') as r:
        for line in r:
            Homo_list.append(line.strip())

    bed_list = os.listdir(bed_dir)
    """ # 不需要物种名，这边重新运行
    with open(selected_species_id, 'r') as r:
        for line in r:
            c = line.strip().split("\t")
            if c[0] in Homo_list:
                gene_sp_dic[c[0]] = c[1]
    """

    # print(Homo_list)
    os.system("mkdir big_fasta")
    for gene_name in Homo_list:
        # gene_name = line.strip()
        # bed_file_path = gene_name + ".bed"
        # if os.path.exists(bed_file_path):
        #     pass
        # else:
        if gene_name+".bed" in bed_list:
            bed_file_path = os.path.join(bed_dir, gene_name + ".bed")
        else:
            for i in bed_list:
                if gene_name in i:
                    bed_file_path = os.path.join(bed_dir, i)
        # bed_file_path = os.path.join("ALL_bed", gene_name + ".bed")

        gene_list_df = pd.read_csv(bed_file_path, header=None, sep='\t')
        gene_list_old = gene_list_df.iloc[:, 4].tolist()
        set_a = set(gene_list_old)
        set_b = set(Homo_list)
        set_gene_list = set_a - set_b
        gene_list = list(set_gene_list)
        # for g in gene_list:
        #     if g in Homo_list:
        #         gene_list.remove(g)
        #     else:
        #         pass
        print(gene_list)
        get_fa_from_id(all_pep, gene_name + '.fa', gene_list)

        # gene_list.to_csv('tmp.gene', index=False, header=False)
        # get_fa_from_id(all_pep, gene_name + '.fa', 'tmp.gene') # 候选成簇的基因蛋白提取

        go_dic={}
        for index, row in gene_list_df.iterrows():
            if row.iloc[4] in Homo_list:
                pass
            else:
                if row.iloc[3] == "+":
                    go_dic[row.iloc[4]] = 1
                else:
                    go_dic[row.iloc[4]] = -1


        # gene_list = gene_list_df.iloc[:, 4].tolist() # 形成列表
        domtblout_file = os.path.join(domout_DB, gene_name + ".domtblout")
        if os.path.exists(domtblout_file):
            pass 
        else:
            with open(all_hmmout,'r') as r, open(domtblout_file, 'w') as w:
                for line in r:
                    if line.startswith("#"):
                        pass
                    else:
                        content = line.strip().split()
                        # print(content[3])
                        if content[3] in gene_list:
                            w.write("\t".join(content) + "\n")
                        else:
                            pass
            # hmm_sqlite3.select_data(domtblout_file, gene_list, gene_sp_dic[gene_name], HMMSCAN_DB)
            # 从数据库中读取这个蛋白的注释结果
            # hmmscan_cmd = "hmmscan --cpu 0 --domtblout {} --cut_tc {} {}.fa".format(domtblout_file, pfam_A_path, gene_name)
            # os.system(hmmscan_cmd)

        # 对这个基因形成的簇读取为字典
        fasta_dict = fasta_parser(gene_name + ".fa") # 形成字典方便调用
        os.system("mv {}.fa big_fasta/".format(gene_name)) # 存到子文件夹下
        # print(fasta_dict)

        DomainList[gene_name] = []
        DomainNumberGene[gene_name] = []
        CoreGeneIndex[gene_name] = []
        GeneOrientation[gene_name] = []
        BGCs[gene_name] = {}

        with open(domtblout_file, 'r') as r2:
            old_g = ''
            gene_index = 0
            for annotation in r2: # 按行读取
                if annotation.startswith("#"):
                    pass
                else:
                    content = annotation.strip().split("\t")
                    DomainList[gene_name].append(content[1])
                    if content[3] == old_g:
                        domainnumber += 1
                    else:
                        gene_index += 1
                        GeneOrientation[gene_name].append(go_dic[content[3]])

                        domainnumber = 1
                        if old_g == '':
                            pass
                        else:
                            DomainNumberGene[gene_name].append(domainnumber)
                        old_g = content[3]
                        
                    if content[3] == gene_name:
                        if len(CoreGeneIndex[gene_name]) == 1:
                            pass
                        else:
                            CoreGeneIndex[gene_name].append(gene_index)
                    else:
                        pass

                    if content[1] in BGCs[gene_name]:
                        pass 
                    else:
                        BGCs[gene_name][content[1]] = []
                    # BGCs[gene_name][content[0]].append(content[3])
                    BGCs[gene_name][content[1]].append(content[3] + ":"+content[19] + ":" + content[20])

                    domain = content[1]
                    if domain in domain_list:
                        pass 
                    else:
                        domain_list.append(domain)
                    header = content[3]
                    seq = fasta_dict[header]


                    domain_seq_file = open(os.path.join(domains_seq_folder, domain + ".fasta"), 'a') # append 'a'
                    domain_seq_file.write(">{}:{}:{}\n{}\n".format(header, content[19], content[20], seq[int(content[19])-1:int(content[20])]))
                    domain_seq_file.close()
            DomainNumberGene[gene_name].append(domainnumber)
    # 对所有domain进行alignment
    for domain in domain_list: # 生产stk文件
        domain_file = os.path.join(domains_seq_folder, domain + ".fasta")
        hmmfetch_pars = ["hmmfetch", pfam_A_path + ".h3m", domain]
        proc_hmmfetch = subprocess.Popen(hmmfetch_pars, stdout=subprocess.PIPE, shell=False)
    
        domain_file_stk = domain_file[:-6]+".stk"
        hmmalign_pars = ["hmmalign", "-o", domain_file_stk, "-", domain_file]
        proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout, stdout=subprocess.PIPE, shell=False)
        
        proc_hmmfetch.stdout.close()
        proc_hmmalign.communicate()[0]
        proc_hmmfetch.wait()

        # print(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))

        stockholm_parser(domain_file_stk) # 从stk生产algn

    AlignedDomainSequences = {}

    for domain in domain_list:
        aligned_file = os.path.join(domains_seq_folder, domain + ".algn")
        # with open(aligned_file, "r") as aligned_file_handle:
        fasta_dict = fasta_parser(aligned_file)
        for header in fasta_dict:
            if "|" in header:
                AlignedDomainSequences[header.split("|")[1]] =fasta_dict[header]
            else:
                AlignedDomainSequences[header] = fasta_dict[header]

    # 循环遍历值列表
    Jaccard_Index = {}
    DSS = {}
    rdssna = {}
    rdss = {}
    for gene_a in Homo_list:
        Jaccard_Index[gene_a] ={}
        DSS[gene_a] = {}
        rdssna[gene_a] = {}
        rdss[gene_a] = {}
        for gene_b in Homo_list:
            domain_list_A = DomainList[gene_a]
            domain_list_B = DomainList[gene_b]
            dcg_a = DomainNumberGene[gene_a]
            dcg_b = DomainNumberGene[gene_b]

            core_pos_a = CoreGeneIndex[gene_a]
            core_pos_b = CoreGeneIndex[gene_b]
            
            # go = "gene orientation"
            go_a = GeneOrientation[gene_a]
            go_b = GeneOrientation[gene_b]
            
            dist, jaccard, dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB, seedLength, reverse = cluster_distance_lcs(gene_a, gene_b, domain_list_A,
                domain_list_B, dcg_a, dcg_b, core_pos_a, core_pos_b, go_a, go_b, 'xxx', BGCs, AlignedDomainSequences)

            Jaccard_Index[gene_a][gene_b] = jaccard
            DSS[gene_a][gene_b] = dss
            rdssna[gene_a][gene_b] = rdssna
            rdss[gene_a][gene_b] = rDSS
        # Homo_list.remove(gene_a)
    Jaccard_df = pd.DataFrame(Jaccard_Index)
    DSS_df = pd.DataFrame(DSS)
    Jaccard_df.to_csv("Jaccard_Index.txt", sep ='\t', header=True, index=True, index_label='Index')
    DSS_df.to_csv("DSS.txt", sep='\t', header=True, index=True, index_label='Index')
    # rDSSna_df = pd.DataFrame(rdssna)
    # rDSS_df = pd.DataFrame(rdss)
    # rDSSna_df.to_csv("rDSSna.txt", sep ='\t', header=True, index=True)
    # rDSS_df.to_csv("rDSS.txt", sep='\t', header=True, index=True)
    merge_files("Jaccard_Index.txt", "DSS.txt", "Index_hebing.txt")

DSS_and_JaccardIndex(homo_list, bed_dir, all_pep, all_hmmout, pfam_A_path)