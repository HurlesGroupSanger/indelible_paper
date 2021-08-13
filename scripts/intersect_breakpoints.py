#!/usr/bin/env python3

from intervaltree import Interval, IntervalTree
import sys

bait_file = open("rawdata/benchmark/SureSelect_baits_padded.bed", 'r')

baits = {}

for bait in bait_file:
    current_rec = bait.rstrip()
    current_rec = current_rec.split("\t")

    current_chr = current_rec[0]
    start = int(current_rec[1])
    end = int(current_rec[2])
    if current_chr in baits:
        baits.get(current_chr).addi(start, end, '1')
    else:
        current_tree = IntervalTree()
        current_tree.addi(start, end, '1')
        baits[current_chr] = current_tree


bed_file = open(sys.argv[1], 'r')

for rec in bed_file:

    current_rec = rec.rstrip()
    data = current_rec.split("\t")
    chrom = data[0]
    start = int(data[1])
    end = int(data[2])
    ref = len(data[3])
    alt = len(data[4])

    var_len = alt - ref

    if chrom in baits:
        start_overlap = baits[chrom].overlaps(start, start + 1)
        end_overlap = baits[chrom].overlaps(end, end + 1)

        if start_overlap is True or end_overlap is True:

            print("%s\t%i\t%i\t%i" % (chrom,start,end,var_len))
        
