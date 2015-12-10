#!/usr/bin/env python3

import os
import glob
import argparse
import statistics
from statistics import median

class Window:
    def __init__(self, scaffold, start, end, gc, readdepth):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.gc = int(gc)
        self.readdepth = float(readdepth)

    def __repr__(self):
        output = "{}\t{}\t{}\t{}\t{}".format(self.scaffold,self.start,self.end,self.gc, round(self.readdepth))
        return(output)

    def adjustreaddepth(self, f):
        self.adjustedreaddepth = f[self.gc]*self.readdepth

parser=argparse.ArgumentParser(description='''Calculate read depths adjusted for GC, genome-wide median read depth, and adjusted genome size
    -w windows
    -o output
    -s size of windows''')

parser.add_argument('-w', '--windows', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-s', '--size', type=int, required=True)

args = parser.parse_args()

print("Loading windows...")
windows = []
gc_median_readdepths = {}
with open(args.windows) as w:
    for window in w:
        scaffold, start, end, gc, readdepth = window.rstrip().split("\t")
        windows.append(Window(scaffold, start, end, gc, readdepth))
        if windows[-1].gc not in gc_median_readdepths:
            gc_median_readdepths[windows[-1].gc] = []
        gc_median_readdepths[windows[-1].gc].append(windows[-1].readdepth)

gc_medians = {}
for gc in gc_median_readdepths:
    gc_medians[gc] = (len(gc_median_readdepths[gc]), median(gc_median_readdepths[gc]))

genome_median = median([w.readdepth for w in windows])
print("Genome median readdepth: {}".format(genome_median))

f = {}
for gc in sorted(gc_medians):
    f[gc] = genome_median / gc_medians[gc][1]
    print("{}\t{}\t{}\t{:.2f}".format(gc, gc_medians[gc][0], gc_medians[gc][1], f[gc]))

windowbases = 0
adjustedbases = 0
unadjustedbases = 0
with open(args.output, 'w') as o:
    o.write("Scaffold\tStart\tEnd\tGC\tReadDepth\tAdjustedReadDepth\n")
    for window in windows:
        window.adjustreaddepth(f)
        o.write("{}\t{}\n".format(repr(window), round(window.adjustedreaddepth)))
        windowbases += args.size
        unadjustedbases += window.readdepth / genome_median * args.size
        adjustedbases += window.adjustedreaddepth / genome_median * args.size

print("Genome size from covered windows: {}".format(windowbases))
print("Genome size by read depth, unadjusted for GC: {}".format(int(unadjustedbases)))
print("Genome size by read depth, adjusted for GC: {}".format(int(adjustedbases)))
