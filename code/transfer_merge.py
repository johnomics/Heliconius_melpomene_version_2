#!/usr/bin/python3 -u

import os
import sys
import re
import pathlib
import argparse
import sqlite3 as sql
from collections import defaultdict, namedtuple
from operator import itemgetter

Haplotype = namedtuple('Haplotype', 'name start end strand trail')

class Part():
    def __init__(self, oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, haplist):
        self.oldname = oldname
        self.oldstart = int(oldstart)
        self.oldend = int(oldend)
        self.newname = newname
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = int(strand)
        self.parttype = parttype
        self.haplotype = None
        if haplist:
            self.haplotype = Haplotype(haplist[0], haplist[1], haplist[2], haplist[3], haplist[4])

    @property
    def length(self):
        return self.oldend - self.oldstart + 1
    
    def __repr__(self):
        out = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.oldname, self.oldstart, self.oldend, self.newname, self.newstart, self.newend, self.strand, self.parttype)
        if self.haplotype:
            out += '\t{}\t{}\t{}\t{}\t{}'.format(self.haplotype.name, self.haplotype.start, self.haplotype.end, self.haplotype.strand, self.haplotype.trail)
        return out


def load_genome(mergedgenome):
    genome = defaultdict(list)
    try:
        with open(mergedgenome, 'r') as g:
            for line in g:
                oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, *args = line.rstrip().split('\t')
                if strand == '0':
                    strand = '1'
                genome[oldname].append(Part(oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, args))
    except IOError:
        print("Failed to load genome file {}".format(mergedgenome))
        sys.exit()
    
    return genome

def get_parts(scaffold, start, end, genome):

    parts = []
    for part in genome[scaffold]:
        if part.oldstart > end or part.oldend < start:
            continue

        slice_start = max(part.oldstart, start)
        slice_end = min(part.oldend, end)
        slice_length = slice_end - slice_start + 1

        comment = '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(part.oldname, part.oldstart, part.oldend, start, end, slice_start, slice_end)
        if part.parttype in ['haplotype', 'gap']:
            newstart = part.newstart
            newend = part.newend
        else:
            start_offset = slice_start - part.oldstart
            end_offset = part.oldend - slice_end
            if part.strand == 1:
                newstart = part.newstart + start_offset
                newend = part.newend - end_offset
            elif part.strand == -1:
                newstart = part.newstart + end_offset
                newend = part.newend - start_offset
        newlength = newend - newstart + 1
        parts.append([part.newname, newstart, newend, newlength, part.strand, part.parttype, comment, part.haplotype])
    return parts

def order(start, end):
    if start < end:
        return start, end
    else:
        return end, start

def get_trailhap(name, start, end, new):
    happarts = get_parts(name, start, end, new)
    trail = []
    for p in happarts:
        if 'removed' in p[5] and p[7] is not None:
            trail.append('{}:{}:{}'.format(p[7].name, p[7].start, p[7].end))
        else:
            trail.append('{}:{}:{}'.format(p[0],p[1],p[2]))
    return ','.join(trail)

def transfer_genome(new, draft, output, prefix):

    transfer = []
    for scaffold in sorted(draft):
        for part in draft[scaffold]:
            newstart, newend = order(part.newstart, part.newend)

            if part.parttype in ['haplotype','gap']:
                if part.haplotype.trail == '-':
                    trailhap = get_trailhap(part.newname, newstart, newend, new)
                else:
                    trailhaps = []
                    for trailhap in part.haplotype.trail.split(','):
                        name, start, end = trailhap.split(':')
                        trailhap = get_trailhap(name, int(start), int(end), new)
                        if trailhap:
                            trailhaps.append(trailhap)
                    trailhap = ','.join(trailhaps)
                part.haplotype = Haplotype(part.haplotype.name, part.haplotype.start, part.haplotype.end, part.haplotype.strand, trailhap)
                transfer.append(part)
                continue

            parts = get_parts(part.newname, newstart, newend, new)
            if not parts:
                transfer.append(part)
                continue

            start = part.oldstart if part.strand == 1 else part.oldend
            done = {}
            for p in parts:
                if (p[0],p[1]) in done:
                    continue

                origname, origstart, origend, scfstart, scfend, slicestart, sliceend = p[6].split('\t')
                slicestart, sliceend = int(slicestart), int(sliceend)
                offset = sliceend - slicestart
                end = start + offset if part.strand == 1 else start - offset
                outstart, outend = order(start, end)
                transfer.append(Part(part.oldname, outstart, outend, p[0], p[1], p[2], part.strand * p[4], p[5], p[7]))

                # Deal with gaps that appear within haplotypes
                inner_parts = []
                for ip in parts:
                    if ip == p:
                        continue
                    if p[0] == ip[0] and p[1] <= ip[1] and ip[2] <= p[2]:
                        i_origname, i_origstart, i_origend, i_scfstart, i_scfend, i_slicestart, i_sliceend = ip[6].split('\t')
                        i_slicestart, i_sliceend = int(i_slicestart), int(i_sliceend)
                        i_offset = i_sliceend - i_slicestart
                        i_start = start + (i_slicestart - slicestart)
                        i_end   = i_start + i_offset if part.strand == 1 else i_start - i_offset
                        i_outstart, i_outend = order(i_start, i_end)
                        transfer.append(Part(part.oldname, i_outstart, i_outend, ip[0], ip[1], ip[2], part.strand * ip[4], ip[5], ip[7]))
                        done[(ip[0],ip[1])] = True
                        
                start = end + part.strand
                done[(p[0],p[1])] = True

    for scaffold in new:
        if prefix and scaffold.startswith(prefix):
            for part in new[scaffold]:
                transfer.append(part)
    
    try:
        with open(output, 'w') as o:
            transfer.sort(key=lambda x: (x.oldname, x.oldstart))
            for part in transfer:
                o.write(repr(part) + "\n")
    except IOError:
        print("Can't open output file", output)
        sys.exit()


def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -n new_genome
        -d draft
        -o output
        -p prefix
        ''')

    parser.add_argument('-n', '--new_genome', type=str, required=True)
    parser.add_argument('-d', '--draft', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-p', '--prefix', type=str, required=False)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    new = load_genome(args.new_genome)
    draft = load_genome(args.draft)
    
    transfer_genome(new, draft, args.output, args.prefix)