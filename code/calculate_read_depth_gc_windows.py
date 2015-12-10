#!/usr/bin/env python3

import os
import glob
import argparse
from Bio import SeqIO
from math import floor
from Bio.SeqUtils import GC
from statistics import median
from collections import defaultdict

class Window:
    def __init__(self, scaffold, start, end):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.gc = int(GC(sequences[self.scaffold][self.start:self.end].seq))
        self.readdepths = defaultdict(int)

    def __repr__(self):
        output = "{}\t{}\t{}\t{}".format(self.scaffold,self.start,self.end,self.gc)
        if hasattr(self, 'median'):
            output += "\t{}".format(self.median)
        return(output)
    
    def addreaddepth(self, readdepth):
        self.readdepths[readdepth] += 1
    
    def calculate_median_readdepth(self):
        medianlists = [[rd] * self.readdepths[rd] for rd in self.readdepths]
        
        self.median = median([item for sublist in medianlists for item in sublist]) if medianlists else 0

parser=argparse.ArgumentParser(description='''Calculate read depth and GC content from bedtools perbase read depth and windows in BED format
    -b windows to test
    -c read depth per base
    -f reference fasta
    -o output prefix
    -w windowsize
    -s scaffold prefix''')

parser.add_argument('-b', '--bed', type=str, required=True)
parser.add_argument('-c', '--readdepth', type=str, required=True)
parser.add_argument('-f', '--fasta', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-w', '--windowsize', type=int, required=True)
parser.add_argument('-s', '--scaffoldprefix', type=str, default=None)

args = parser.parse_args()

print("Loading sequences...")
sequences = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

print("Loading windows...")
windows = defaultdict(lambda:defaultdict(Window))
with open(args.bed) as b:
    for window in b:
        scaffold, start, end = window.rstrip().split("\t")
        if args.scaffoldprefix and not scaffold.startswith(args.scaffoldprefix):
            continue
        start, end = int(start), int(end)
        if end - start != args.windowsize:
            continue
        windows[scaffold][start] = Window(scaffold, start, end)

print("Loading read depths...")
win = 0
curscaffold = ''
scaffolds_done = False
with open(args.readdepth) as r:
    for base in r:
        scaffold, position, readdepth = base.rstrip().split("\t")
        if args.scaffoldprefix and not scaffold.startswith(args.scaffoldprefix):
            if scaffolds_done:
                break
            else:
                continue
        scaffolds_done = True
        if scaffold != curscaffold or curscaffold == '':
            print(scaffold)
            curscaffold = scaffold
        position, readdepth = int(position)-1, int(readdepth)
        
        if scaffold in windows:
            start = floor(position / args.windowsize) * args.windowsize
            if start in windows[scaffold]:
                windows[scaffold][start].addreaddepth(readdepth)

with open(args.output + ".windows.out", 'w') as o:
    for scaffold in sorted(windows):
        for start in sorted(windows[scaffold]):
            windows[scaffold][start].calculate_median_readdepth()
            o.write(repr(windows[scaffold][start]) + "\n")
