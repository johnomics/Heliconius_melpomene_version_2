#!/usr/bin/python3 -u

import os
import sys
import argparse
import gzip
from collections import defaultdict
from copy import deepcopy

class Overlap():
    def __init__(self, name, size, start, end, length, strand = 0, my_id = 0, new_name=''):
        self.name = name
        self.size = int(size)
        self.start = int(start)+1
        self.end = int(end)
        self.length = int(length)
        self.strand = int(strand)
        self.id = int(my_id)
        self.new_name = new_name
    
    def __repr__(self):
        return '{}:{}-{} {}bp {} {}'.format(self.name, self.start, self.end, self.length, self.strand, self.new_name)

class Part():
    def __init__(self, line, num=-1, prefix = ""):
        self.scaffold_id, self.scaffold_len, self.sub_scaffold_id, self.new_portion_id, self.active_portion, \
        old_scaffold1_name, old_scaffold1_size, old_scaffold1_id, start1, end1, strand1, len1, \
        old_scaffold2_name, old_scaffold2_size, old_scaffold2_id, start2, end2, strand2, len2, \
        self.connection, self.score, self.sc_Ns, self.tsc_LCs, self.qsc_Ns, self.qsc_LCs, \
        self.active_portion_updated, self.connection_updated, self.active_portion_manual, self.connection_manual \
            = line.rstrip().split('\t')

        self.connection = int(self.connection)
        
        if num >= 0:
            self.scaffold_name = prefix + 'Sc{:07d}'.format(num)
        self.scaffold1 = Overlap(old_scaffold1_name, old_scaffold1_size, start1, end1, len1, strand1, old_scaffold1_id)
        self.scaffold2 = Overlap(old_scaffold2_name, old_scaffold2_size, start2, end2, len2, strand2, old_scaffold2_id)
        
    def __repr__(self):
        return '{}\t{}'.format(self.scaffold_name, self.scaffold_len)

class OutPart:
    def __init__(self, name, start, end, new_name, new_start, new_end, strand, connection, parttype, haplotype=None):
        self.name = name
        self.start = start
        self.end = end
        self.new_name = new_name
        self.new_start = new_start
        self.new_end = new_end
        self.strand = strand
        self.connection = connection
        self.type = parttype
        self.haplotype = haplotype
    
    @property
    def length(self):
        return self.end - self.start + 1

    def __repr__(self):
        out = '{}\t{}\t{}\t{}\t{}\t{}'.format(self.end, self.new_name, self.new_start, self.new_end, self.strand, self.type)
        if self.haplotype:
            out += '\t' + self.haplotype.haplorepr()
        return out
    
    def haplorepr(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.new_name, self.new_start, self.new_end, self.strand, '-')

def load_scaffolds(prefix, old_genome, new_genome):

    new_scaffolds_file = "hm.new_scaffolds_updated"
    if not os.path.isfile(new_scaffolds_file):
        new_scaffolds_file = "hm.new_scaffolds"
        if not os.path.isfile(new_scaffolds_file):
            print("Can't find " + new_scaffolds_file)
            sys.exit()
    
    try:
        with open(new_scaffolds_file) as s:
            num = 0
            curid = ''
            curstart = 1
            for line in s:
                if line.startswith('#'):
                    continue

                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                scaffold_id = f[0]
                if not curid:
                    curid = scaffold_id
                
                if scaffold_id != curid:
                    curid = scaffold_id
                    curstart = 1
                    num += 1

                part = Part(line, num, prefix)

                active_portion = part.active_portion_manual if part.active_portion_manual != '0' else part.active_portion
                scaffold = part.scaffold1 if active_portion == '1' else part.scaffold2
                curend = curstart + scaffold.length - 1

                new_genome[part.scaffold_name][curstart] = OutPart(part.scaffold_name, curstart, curend, scaffold.name, scaffold.start, scaffold.end, scaffold.strand, part.connection, 'merged')
                old_genome[scaffold.name][scaffold.start] = OutPart(scaffold.name, scaffold.start, scaffold.end, part.scaffold_name, curstart, curend, scaffold.strand, part.connection, "active")
                hapscaffold = part.scaffold1 if active_portion == '2' else part.scaffold2
                if hapscaffold.name != "0":
                    old_genome[hapscaffold.name][hapscaffold.start] = OutPart(hapscaffold.name, hapscaffold.start, hapscaffold.end, part.scaffold_name, curstart, curend, 
                                                                              hapscaffold.strand, part.connection, "haplotype",
                                                                              new_genome[part.scaffold_name][curstart])

                curstart = curend + 1

    except IOError:
        print("Failed to load scaffolds file {}".format(scaffolds_path))
        sys.exit()
    
    new_genome = reorder_genomes(old_genome, new_genome, prefix)
    return new_genome

def load_refined():
    
    refined_types = defaultdict(str)
    refined_file = "../_hm.unincorpRefiner.log"
    if not os.path.isfile(refined_file):
        print("Can't find _hm.unincorpRefiner.log, no refined types will be output")
        return refined_types

    try:
        with open(refined_file) as r:
            for line in r:
                if not (line.startswith('xp') or line.startswith('xf')):
                    continue
                
                new_name, old_name, size, start, length, fulltype, ns, lowcase, ali_len, ali_cov, retained = line.rstrip().split('\t')
                if retained == "1":
                    refined_types[new_name] = "retained"
                else:
                    refined_types[new_name] = "removed"

            return refined_types

    except IOError:
        print("Failed to load refined file {}".format(refined_file))
        sys.exit()

def load_refined_hits(refined_types, prefix):
    refined_hits = defaultdict(lambda: defaultdict(list))
    refined_file = "unpaired.tbest.net.gz"
    if not os.path.isfile(refined_file):
        print("Can't find unpaired.tbest.net.gz, no refined hits will be output")
        return refined_hits

    try:
        with gzip.open(refined_file, 'rt') as r:
            red_name = ""
            for line in r:
                if not (line.startswith('net') or line.startswith(' fill')):
                    continue

                if line.startswith('net'):
                    net, red_name, length = line.rstrip().split(' ')
                    red_name = prefix + red_name
                elif line.startswith(' fill'):
                    x, fill, red_start, red_length, new_name, strand, new_start, new_length, *args = line.rstrip().split(' ')
                    red_start, red_length, new_start, new_length = int(red_start), int(red_length), int(new_start), int(new_length)
                    red_end = red_start + red_length
                    new_end = new_start + new_length
                    red_start = red_start + 1
                    new_start = new_start + 1
                    new_name = prefix + new_name
                    part = OutPart(red_name, red_start, red_end, new_name, new_start, new_end, strand, 0, refined_types[red_name])
                    refined_hits[red_name][new_name].append(part)
    except IOError:
        print("Failed to load refined hits file {}".format(refined_file))
    
    for red_name in refined_hits:
        for new_name in refined_hits[red_name]:
            if len(refined_hits[red_name][new_name]) > 1:
                min_red_start = min_new_start = 100000000
                max_red_end = max_new_end = 0
                strand = 0
                parttype = ''
                for part in refined_hits[red_name][new_name]:
                    strand = part.strand
                    parttype = part.type
                    if part.start < min_red_start:
                        min_red_start = part.start
                    if part.end > max_red_end:
                        max_red_end = part.end
                    if part.new_start < min_new_start:
                        min_new_start = part.new_start
                    if part.new_end > max_new_end:
                        max_new_end = part.new_end

                refined_hits[red_name][new_name] = [OutPart(red_name, min_red_start, max_red_end, new_name, min_new_start, max_new_end, strand, 0, parttype)]

    return refined_hits

def load_unpaired(prefix, old_genome, new_genome, lengths):
    
    refined_types = load_refined()
    refined_hits = load_refined_hits(refined_types, prefix)
    
    if not os.path.isfile("hm.unpaired_updated"):
        print("Can't find hm.unpaired_updated")
        sys.exit()
    
    try:
        with open("hm.unpaired_updated") as u:
            for line in u:
                if line.startswith('#'):
                    continue
                
                f = line.rstrip().split('\t')
                if len(f) <= 1:
                    continue
                new_name, old_name, size, start, length, fulltype, ns, lcs = f
                size, start, length = int(size), int(start), int(length)

                refined_type = fulltype
                if new_name in refined_types:
                    refined_type = refined_types[new_name]

                new_name = prefix + new_name
                
                lengths[new_name] = (length, length)

                overlap = Overlap(old_name, size, start, start+length, length, new_name=new_name)
                new_start = 1
                new_end = overlap.length
                
                new_genome[new_name][new_start] = OutPart(new_name, new_start, new_end, overlap.name, overlap.start, overlap.end, overlap.strand, 0, refined_type)
                old_genome[overlap.name][overlap.start] = OutPart(overlap.name, overlap.start, overlap.end, new_name, new_start, new_end, overlap.strand, 0, refined_type)
                if new_name in refined_hits and refined_type == 'removed':
                    ordered_hits = sorted(refined_hits[new_name], key=lambda x: refined_hits[new_name][x][0].length, reverse=True)
                    old_genome[overlap.name][overlap.start].haplotype = refined_hits[new_name][ordered_hits[0]][0]
                
    except IOError:
        print("Failed to load unpaired file {}".format(unpaired_path))

class Node:
    def __init__(self, target_name, query_name, node_id, tstart, tend, qstart, qend, strand):
        self.tname = target_name
        self.qname = query_name
        self.node_id = node_id
        self.tstart = tstart
        self.tend = tend
        self.qstart = qstart
        self.qend = qend
        self.strand = strand

    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.node_id, self.tname, self.tstart, self.tend, self.qname, self.qstart, self.qend, self.strand)
        
def load_nodes(scaffold_prefix):
    nodes = defaultdict(Node)
    nodes_filename = "hm.nodes_updated"
    if not os.path.isfile(nodes_filename):
        nodes_filename = "hm.nodes"
        if not os.path.isfile(nodes_filename):
            print("Can't find " + nodes_filename)
            sys.exit()

    try:
        with open(nodes_filename, 'r') as nodes_file:
            for line in nodes_file:
                if line.startswith('#'):
                    continue
                f = line.rstrip().split('\t')
                if len(f) == 1:
                    continue
                target_name = f[0]
                if target_name.startswith(scaffold_prefix):
                    continue
                query_name, node_id, tstart, tend, qstart, qend, strand = f[1], int(f[2]), int(f[5])+1, int(f[7]), int(f[10])+1, int(f[12]), int(f[9])
                nodes[node_id] = Node(target_name, query_name, node_id, tstart, tend, qstart, qend, strand)
    except IOError:
        print("Failed to load nodes file")

    return nodes


class Gap:
    def __init__(self, scaffold, ali_len, gap_start, gap_end, new_name, node_id, strand, old_start, old_end, new_start, new_end):
        self.scaffold  = scaffold
        self.start  = gap_start
        self.end    = gap_end
        self.length = self.end - self.start + 1
        self.new_name = new_name
        self.node_id  = node_id
        self.strand = strand
        self.old_start = old_start
        self.old_end   = old_end
        self.new_start = new_start
        self.new_end   = new_end
        self.old_length = self.old_end - self.old_start + 1
        self.new_length = new_end - new_start + 1
        self.fixed = False
        
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.strand, self.length, 
                                                               self.old_start, self.old_end, self.old_length,
                                                               self.new_name, self.new_start, self.new_end, self.new_length)

def load_gaps(nodes):
    gaps = defaultdict(lambda: defaultdict(Gap))
    ngap_file = "hm.n_gaps"
    if not os.path.isfile(ngap_file):
        return gaps

    try:
        with open(ngap_file, 'r') as gapfile:
            for line in gapfile:
                if line.startswith('#'):
                    continue
                f = line.rstrip().split('\t')
                alignment_length = int(f[8])
                gap_pc_remaining = int(f[21])
                if alignment_length < 4000 or gap_pc_remaining > 10:
                    continue
                scaffold, gap_start, gap_end, node_id, strand, old_start, old_end, new_start, new_end = f[0], int(f[3])+1, int(f[4]), int(f[6]), int(f[10]), int(f[17])+1, int(f[18]), int(f[19])+1, int(f[20])
                fill_scaffold = nodes[node_id].qname
                gaps[scaffold][old_start] = Gap(f[0], alignment_length, gap_start, gap_end, fill_scaffold, node_id, strand, old_start, old_end, new_start, new_end)

    except IOError:
        print("Failed to load gap file")

    return gaps

def load_merge(prefix, mergedir, old_genome, new_genome):

    results = mergedir + "/genome.genomex.result"
    if not os.path.isdir(results):
        print("Results argument is not a directory!")
        sys.exit()
    
    try:
        os.chdir(results)
    except OSError:
        print("Could not change directory to results directory!")
        sys.exit()

    lengths = load_lengths()

    new_genome = load_scaffolds(prefix, old_genome, new_genome)
    
    load_unpaired(prefix, old_genome, new_genome, lengths)
    
    return new_genome, lengths

def collapse(genome):
    for scaffold in genome:
        starts = sorted(genome[scaffold].keys())
        for start_i in starts:
            if start_i not in genome[scaffold]:
                continue
            for start_j in starts:
                if start_j not in genome[scaffold] or start_j <= start_i:
                    continue
                g_i = genome[scaffold][start_i]
                g_j = genome[scaffold][start_j]
                if g_i.new_name == g_j.new_name and g_i.strand == g_j.strand and g_i.type == g_j.type:
                    if g_i.end + 1 == start_j:
                        if g_i.strand == 1 and g_i.new_end + 1 == g_j.new_start:
                            g_i.end = g_j.end
                            g_i.new_end = g_j.new_end
                            del genome[scaffold][start_j]
                        elif g_i.strand == -1 and g_i.new_start == g_j.new_end + 1:
                            g_i.end = g_j.end
                            g_i.new_start = g_j.new_start
                            del genome[scaffold][start_j]

def reorder_genomes(old, new, prefix):
    lengths = {}
    for scaffold in new:
        length = sum([new[scaffold][start].length for start in new[scaffold]])
        lengths[scaffold] = length
    
    newnames = {}
    num = 0
    new_new_genome = defaultdict(lambda: defaultdict(OutPart))
    for scaffold in sorted(lengths, key = lambda x: lengths[x], reverse=True):
        newname = prefix + 'Sc{:07d}'.format(num)
        newnames[scaffold] = newname
        new_new_genome[newname] = new[scaffold]
        num += 1

    for old_scaffold in old:
        for start in old[old_scaffold]:
            old[old_scaffold][start].new_name = newnames[old[old_scaffold][start].new_name]

    return new_new_genome

def update_gap(gap, scaffold, partstart, old, new, lengths):

    new_scaffold = old[scaffold][partstart].new_name

    if (   lengths[new_scaffold][0] == lengths[new_scaffold][1]         # No change in length
        or scaffold == gap.new_name                                     # Skip gaps filled with different parts of the same scaffold, as HaploMerger also skips these
        or old[scaffold][partstart].type not in ['active', 'retained']  # Only fill gaps in regions contained in the genome, ignore haplotypes
        or old[scaffold][partstart].connection < 0):                    # A small number of gaps have connection = -2, no connection, so no gap can be filled; a node is present but rejected before hm.new_scaffolds
        gap.fixed = True
        return

    part = deepcopy(old[scaffold][partstart])

    offset = gap.new_length - gap.old_length
    if part.strand in [1,0]:    # Forward
        new_gap_start = part.new_start + gap.old_start - partstart
        new_gap_end   = new_gap_start + gap.new_length - 1
        new_end = new_gap_end + part.end - gap.old_end
    elif part.strand == -1:     # Reverse
        new_gap_start = part.new_end - gap.old_end + part.start
        new_gap_end   = part.new_end - gap.old_start + part.start + offset
        new_end       = part.new_end + offset

    for oldscaffold in old:
        for oldstart in old[oldscaffold]:
            if oldscaffold == scaffold and oldstart == partstart:
                continue
            oldpart = old[oldscaffold][oldstart]
            if oldpart.new_name == part.new_name:
                if oldpart.new_start >= new_gap_start:
                    oldpart.new_start += offset
                if oldpart.new_end >= new_gap_start:
                    oldpart.new_end += offset
            
    del(old[scaffold][partstart])

    # Output pre-gap section
    if part.strand in [1,0]:
        old_start, old_end = partstart, gap.old_start-1
    elif part.strand == -1:
        old_start, old_end = gap.old_end+1, part.end
    old[scaffold][old_start] = OutPart(scaffold, old_start, old_end, part.new_name, part.new_start, new_gap_start - 1, part.strand, part.connection, part.type)

    # Output gap and removed section
    new_gap_strand = int(gap.strand) * int(part.strand)
    old[gap.new_name][gap.new_start] = OutPart(gap.new_name, gap.new_start, gap.new_end, part.new_name, new_gap_start, new_gap_end, new_gap_strand, part.connection, 'active')

    old[scaffold][gap.old_start] = OutPart(scaffold, gap.old_start, gap.old_end, part.new_name, new_gap_start, new_gap_end, part.strand, part.connection, 'gap',
                                            OutPart(scaffold, gap.old_start, gap.old_end, gap.new_name, gap.new_start, gap.new_end, gap.strand, part.connection, 'gap'))

    # Output post-gap section
    if part.strand == 1:
        old_start, old_end = gap.old_end+1, part.end
    elif part.strand == -1:
        old_start, old_end = partstart, gap.old_start-1
    old[scaffold][old_start] = OutPart(scaffold, old_start, old_end, part.new_name, new_gap_end+1, new_end, part.strand, part.connection, part.type)

    gap.fixed = True


def fill_gaps(old, new, lengths, args):
    
    nodes = load_nodes(args.scaffolding)
    gaps = load_gaps(nodes)

    for scaffold in sorted(old):
        if scaffold.startswith(args.scaffolding):
            continue
        i = 0
        while (i < len(old[scaffold])):
            recheck = False
            starts = sorted(old[scaffold])
            start = starts[i]
            end = old[scaffold][start].end
            if scaffold in gaps:
                for gap_start in sorted(gaps[scaffold]):
                    if gaps[scaffold][gap_start].fixed:
                        continue
                    gap_end = gaps[scaffold][gap_start].end
                    if start <= gap_start <= end or start <= gap_end <= end:
                        update_gap(gaps[scaffold][gap_start], scaffold, start, old, new, lengths)
                        recheck = True
                        break
            if not recheck:
                i += 1

def load_lengths():
    lengths = {}
    with open("../_hm.haploMerger_updated.log") as hm:
        for line in hm:
            if not line.startswith("Sc"):
                continue
            f = line.rstrip().split('\t')
            if len(f) <= 1:
                continue
            scaffold, start_id, end_id, size, refined_size = f
            scaffold = args.prefix+scaffold
            lengths[scaffold] = (int(size), int(refined_size))
    return(lengths)

def sanity_check(old, lengths, args):

    new_lengths = {}
    for scaffold in old:
        for start in old[scaffold]:
            part = old[scaffold][start]
            if part.new_name not in new_lengths:
                new_lengths[part.new_name] = 0
            if part.new_end > new_lengths[part.new_name]:
                new_lengths[part.new_name] = part.new_end

    for scaffold in sorted(lengths):
        if lengths[scaffold][1] != new_lengths[scaffold]:
            print("Length error\t{}\t{}\t{}\t{}".format(scaffold, lengths[scaffold][0], lengths[scaffold][1], new_lengths[scaffold]))


def output_genome(suffix, genome, mergedir):
    
    collapse(genome)

    try:
        file = mergedir + "_" + suffix + ".tsv"
        f = open(file, 'w')
    except IOError:
        print("Can't open " + suffix + " genome output file!")
        sys.exit()

    for scaffold in sorted(genome):
        for start in sorted(genome[scaffold]):
            f.write("{}\t{}\t{}\n".format(scaffold, start, genome[scaffold][start]))

def get_args():
    parser = argparse.ArgumentParser(description='''Output map of new scaffolds to old based on HaploMerger results
        -m HaploMerger directory
        -p output prefix
        -s scaffolding prefix
        ''')

    parser.add_argument('-m', '--mergedir', type=str, required=True)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-s', '--scaffolding', type=str, required=False)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    old_genome = defaultdict(lambda: defaultdict(OutPart))
    new_genome = defaultdict(lambda: defaultdict(OutPart))

    new_genome, lengths = load_merge(args.prefix, args.mergedir, old_genome, new_genome)
    
    if args.scaffolding:
        fill_gaps(old_genome, new_genome, lengths, args)

    sanity_check(old_genome, lengths, args)

    output_genome("old", old_genome, args.mergedir)
