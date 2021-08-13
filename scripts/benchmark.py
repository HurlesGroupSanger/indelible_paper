#!/usr/bin/env python3

from intervaltree import Interval, IntervalTree
import sys

def fill_tree(file_name, vars):

    fh = open(file_name, 'r')
    
    for rec in fh:
    
        current = rec.rstrip()
        data = current.split("\t")

        chrom = data[0]
        start = int(data[1])
        end = int(data[2]) + 1
        length = int(data[3])
        
        var_info = {'length': length,
                    'found': False}
        
        if chrom in vars:
            vars.get(chrom).addi(start, end, var_info)
        else:
            current_tree = IntervalTree()
            current_tree.addi(start, end, var_info)
            vars[chrom] = current_tree

vars = {}

small = fill_tree("rawdata/benchmark/HG002_GRCh37_1_22_v4.1_draft_benchmark.norm.filtered.baits.bed", vars)
large = fill_tree("rawdata/benchmark/HG002_SVs_Tier1_v0.6.filtered.baits.bed", vars)

to_benchmark = sys.argv[1]

fh = open(to_benchmark, 'r')

false_positives = {}

for rec in fh:

    current = rec.rstrip()
    data = current.split("\t")
    chrom = data[0]
    start = int(data[1])
    end = int(data[2])

    current_tree = vars.get(chrom)

    if current_tree is not None:
    
        overlaps = current_tree.overlap(start - 100, end + 100)

        found_overlaps = 0
    
        for interval in overlaps:
            int_start = interval[0]
            int_end = interval[1]
            int_data = interval[2]

            start_dist = min(abs(start - int_start), abs(start - int_end))
            end_dist = min(abs(end - int_start), abs(end - int_end))
            
            if start_dist <= 100 or end_dist <= 100:
                found_overlaps += 1
                int_data['found'] = True
                vars.get(chrom).addi(int_start, int_end, int_data)
            
        if found_overlaps == 0:

            key = "%s_%i" % (chrom, start)
            false_positives[key] = {'chrom': chrom, 'start': start, 'end': end}

for chrom in vars:

    current_tree = vars.get(chrom)

    for obj in current_tree:

        printable = "%s\t%i\t%i\t%i\t%s\tREF" % (chrom, obj[0], obj[1], obj[2]['length'], obj[2]['found'])
        print(printable)
        

for var in false_positives:

    rec = false_positives[var]
    printable = "%s\t%i\t%i\tNA\tNA\tFP" % (rec['chrom'], rec['start'], rec['end'])
    print(printable)
