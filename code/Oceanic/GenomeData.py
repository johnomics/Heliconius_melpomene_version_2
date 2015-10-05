#!/usr/bin/env python3

import sys
import gzip
import re
from os.path import isfile
import sqlite3 as sql
from collections import defaultdict
from termcolor import colored

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from BCBio import GFF

from . import Stats

class GenomeData:
    def __init__(self, args):
        self.haplotypes = {}
        self.refuse = []

        
        self.sequences = self.load_genome(args.fasta)
        
        self.load_annotation(args.gff)
        
        conn = self.open_database(args.database)
        self.db = conn.cursor()
        
        self.blocks = self.load_blocks()
        
        self.errors = self.load_errors(args.errors)

        self.nodes = self.load_nodes(args.nodes)
        
        self.revised, self.revised_fasta, self.revised_tsv, self.revised_orig_tsv, self.revised_db, self.revised_conn = self.open_revised(args.revised)
        
        self.revised_names = {}
        self.revised_count = 1
        
        self.origparts, self.newparts, self.offcuts = self.clean_assembly(args.original, args.prefix)

        self.gapnum = 1
        
        if args.haplomerger:
            print("Loading haplotypes...")
            self.haplotypes = self.load_haplotypes(args.haplomerger)
        else:
            self.haplotypes = {}
        
    def open_revised(self, revised):
        fasta = None
        tsv = None
        orig_tsv = None
        db = None
        rconn = None
        if revised:
            fasta = open(revised+".fa", 'w')
            tsv = open(revised+".tsv", 'w')
            orig_tsv = open(revised+"_orig.tsv", 'w')
            rconn = sql.connect(revised+".db")
            db = rconn.cursor()
            db.execute('drop table if exists scaffold_map')
            db.execute('''create table scaffold_map
                         (chromosome integer, cm real, scaffold text, start integer, end integer, length integer)''')

        return revised, fasta, tsv, orig_tsv, db, rconn

    def open_database(self, dbfile):
        try:
            if isfile(dbfile):
                conn = sql.connect(dbfile)
            else:
                raise IOError
        except IOError:
            print("Can't open database {}!".format(dbfile))
            sys.exit()
    
        return conn

    def load_genome(self, fasta):
        print("Loading genome...")
        try:
            with open(fasta, 'r') as f:
                sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
                for scaffold in sequences:
                    sequences[scaffold].seq = sequences[scaffold].seq.upper()

        except IOError:
            print("Can't load genome from file {}!".format(fasta))
            sys.exit()

        print("Original stats:", Stats.genome([len(sequences[scaffold]) for scaffold in sequences]))

        return sequences


    def load_annotation(self, gff):
        print("Loading annotation...")
        try:
            with open(gff, 'r') as g:
                for feature in g:
                    if 'gene' in feature:
                        scaffold, source, gfftype, start, end, score, strand, phase, attributes = feature.rstrip().split('\t')
                        if scaffold not in self.sequences:
                            continue
                        featurestrand = None
                        if strand == '+':
                            featurestrand = 1
                        elif strand == '-':
                            featurestrand = -1
                        location = FeatureLocation(ExactPosition(start), ExactPosition(end), strand=featurestrand)
                        feature = SeqFeature(location, type=gfftype)
                        self.sequences[scaffold].features.append(feature)
        except IOError:
            print("Can't load annotation from file {}!".format(gff))
            sys.exit()

    def load_blocks(self):
        print("Loading blocks...")
        
        blocks = defaultdict(lambda:defaultdict(Block))
        block_num = 0
        genome_length = 0

        # Load map blocks from database
        for scaffold, start, end, parttype in self.db.execute("select scaffold, start, end, parttype from scaffold_map"):
            if parttype not in ['active', 'retained']:
                continue
            if scaffold not in self.sequences:
                continue
            if end < start:
                continue
            block_num += 1
            genome_length += end - start + 1
            blocks[scaffold][start] = Block(scaffold, start, end)
        
        for scaffold in self.sequences:
            if scaffold not in blocks:
                end = len(self.sequences[scaffold])
                block_num += 1
                genome_length += end
                blocks[scaffold][1] = Block(scaffold, 1, end)

        # Store previous and next blocks
        for scaffold in blocks:
            positions = sorted(blocks[scaffold].keys())
            for i, pos in enumerate(positions):
                if i > 0:
                    blocks[scaffold][pos].prev_block = blocks[scaffold][positions[i-1]].start
                if i < len(positions)-1:
                    blocks[scaffold][pos].next_block = blocks[scaffold][positions[i+1]].start

        print("Loaded {} blocks, length {}".format(block_num, genome_length))
        return blocks

    def load_errors(self, errorfile):
        print("Loading errors...")
        errors = {}
        if not errorfile:
            return errors
    
        try:
            if isfile(errorfile):
                with open(errorfile) as err:
                    for line in err:
                        scaffold, start = line.rstrip().split('\t')
                        errors[(scaffold, int(start))] = -1
            else:
                raise IOError
        except IOError:
            print("Can't open errors file!")
            sys.exit()
    
        return errors

    def load_nodes(self, nodesfile):
        print("Loading nodes...")
                
        nodes = defaultdict(lambda:defaultdict(list))
        if not nodesfile:
            return nodes
        
        try:
            if isfile(nodesfile):
                with open(nodesfile) as n:
                    for line in n:
                        if line.startswith("#"):
                            continue
                        f = line.rstrip().split('\t')
                        if len(f) == 1 or f[0] == f[1]:
                            continue
                        nodes[f[0]][f[1]].append(Node(f[0], f[1], int(f[5])+1, int(f[6]), int(f[7]), int(f[10])+1, int(f[11]), int(f[12]), int(f[9]), self.sequences))
            else:
                raise IOError
        except IOError:
            print("Can't open nodes file", nodesfile)
            sys.exit()


        for tscaffold in list(nodes):
            if tscaffold.endswith("F"):
                continue
            for qscaffold in list(nodes[tscaffold]):
                if qscaffold.endswith("F") and qscaffold in nodes:
                    pbqscaffolds = [x for x in nodes[qscaffold].keys() if not x.endswith("F") and x != tscaffold]
                    if not pbqscaffolds:
                        continue
                    for node in nodes[tscaffold][qscaffold]:
                        for pbqscaffold in pbqscaffolds:
                            if pbqscaffold == tscaffold:
                                continue
                            for pbnode in nodes[qscaffold][pbqscaffold]:
                                nodes[tscaffold][pbqscaffold].append(Node(tscaffold, pbqscaffold, node.tstart, node.tlen, node.tend, pbnode.qstart, pbnode.qlen, pbnode.qend, '0', self.sequences,
                                                                          qscaffold, node.qstart, node.qlen, node.qend, node.direction, pbnode.tstart, pbnode.tlen, pbnode.tend, pbnode.direction))

        for tscaffold in list(nodes):
            if tscaffold.endswith("F"):
                del nodes[tscaffold]
            for qscaffold in list(nodes[tscaffold]):
                if qscaffold.endswith("F"):
                    del nodes[tscaffold][qscaffold]

        return nodes
        
    def clean_assembly(self, tsv, prefix):
        print("Cleaning assembly...")

        origparts = defaultdict(list)
        newparts = defaultdict(list)
        offcuts = defaultdict(lambda: defaultdict(list))
        
        if not tsv:
            return origparts, newparts, offcuts

        try:
            if isfile(tsv):
                with open(tsv) as tsv:
                    for line in tsv:
                        part = OrigPart(line)
                        origparts[part.oldname].append(part)
                        newparts[part.newname].append(part)
        except IOError:
            print("Can't open original TSV file!")
            sys.exit()

        self.delete_scaffolding_only(newparts, origparts, prefix)
        offcuts = self.mark_offcuts(newparts, origparts, prefix)
                
        return origparts, newparts, offcuts
    
    def delete_scaffolding_only(self, newparts, origparts, prefix):
        deleted = 0
        deleted_length = 0
        for new in newparts:
            scaffolding_only = True
            for part in newparts[new]:
                if prefix not in part.oldname:
                    scaffolding_only = False
            if scaffolding_only and new in self.sequences:
                deleted += 1
                deleted_length += len(self.sequences[new])
                del self.sequences[new]
                del self.blocks[new]
                for part in newparts[new]:
                    if part.parttype not in ['haplotype', 'gap']:
                        part.parttype = 'removed'
        
        print("Cleaned up {} scaffolds, length {}".format(deleted, deleted_length))
        
    def mark_offcuts(self, newparts, origparts, prefix):
        offcuts = defaultdict(lambda: defaultdict(list))
        for new in newparts:
            newpart = newparts[new][0]
            if len(newparts[new]) == 1 and prefix not in newpart.oldname and newpart.parttype == 'retained':
                old = newpart.oldname
                if len(origparts[old]) > 1:
                    newpart.parttype = 'offcut'
                    for origpart in origparts[old]:
                        if origpart.newname == new:
                            origpart.parttype = 'offcut'
                        else:
                            if origpart.parttype == 'active':
                                offcuts[new][origpart.newname] = 1
                            elif origpart.parttype == 'haplotype':
                                if origpart.haptrail == '-':
                                    offcuts[new][origpart.newname] = 1
                                else:
                                    for trailhap in origpart.haptrail.split(','):
                                        name, start, end = trailhap.split(':')
                                        offcuts[new][name] = 1

        offcut_length = sum([newparts[x][0].length for x in offcuts])
        print("Found {} offcuts, length {}".format(len(offcuts), offcut_length))
            
        return offcuts

        
    def are_neighbours(self, a, b):
        for ai in a[0], a[1]:
            for bi in b[0], b[1]:
                if abs(int(ai)-int(bi)) == 1:
                    return True
        return False
        
    def add_gap(self, length = 100):
        newgapname = 'Gap' + str(self.gapnum)
        self.blocks[newgapname][1] = Block(newgapname, 1, length)
        self.gapnum += 1
        return self.blocks[newgapname][1]
    
    def load_haplotypes(self, haplomerger):
        haplotypes = {}
        self.load_haplotype(haplotypes, haplomerger, "A")
        self.load_haplotype(haplotypes, haplomerger, "B")
        return haplotypes
    
    def load_haplotype(self, haplotypes, haplomerger, hap):
        hapfile = haplomerger + "/assembly_hap" + hap + ".fa.gz"
        try:
            if isfile(hapfile):
                with gzip.open(hapfile, 'rt') as hf:
                    for line in hf:
                        if line.startswith(">"):
                            f = line[1:].split('_')
                            name = f[0] + '_' + f[1] + '_' + f[2]
                            haplotypes[name] = hap
        except IOError:
            print("Can't open haplotype file!", hapfile)
            sys.exit()

class OrigPart:
    def __init__(self, line):
        self.oldname, self.oldstart, self.oldend, self.newname, self.newstart, self.newend, self.strand, self.parttype, *args = line.rstrip().split('\t')
        self.comment = ''
        self.oldstart, self.oldend, self.newstart, self.newend = int(self.oldstart), int(self.oldend), int(self.newstart), int(self.newend)
        if args:
            self.comment = '\t'.join(args)
            self.hapname, self.hapstart, self.hapend, self.hapstrand, self.haptrail = args
        self.length = self.oldend-self.oldstart+1
    
    def __repr__(self):
        output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.oldname, self.oldstart, self.oldend, self.newname, self.newstart, self.newend, self.strand, self.parttype)
        if self.comment:
            output += "\t{}".format(self.comment)
        return output

class Block:
    def __init__(self, scaffold, start, end, prev_block=0, next_block=0, chromosome=0, cm=-1):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.prev_block = prev_block
        self.next_block = next_block
        self.chromosome = str(chromosome)
        self.cm = cm
    
    @property
    def length(self):
        return self.end - self.start + 1
    
    def __repr__(self):
        return "{}:{}-{} ({}, {}) [{}, {}]".format(self.scaffold, self.start, self.end, self.prev_block, self.next_block, self.chromosome, self.cm)

    def add_marker(self, chromosome, cm):
        self.chromosome = str(chromosome)
        self.cm = cm

class Node:
    def __init__(self, tscaffold, qscaffold, tstart, tlen, tend, qstart, qlen, qend, direction, sequences, *args):
        self.tscaffold = tscaffold
        self.qscaffold = qscaffold
        self.tstart = tstart
        self.tlen = tlen
        self.tend = tend
        self.qstart = qstart
        self.qlen = qlen
        self.qend = qend
        self.direction = int(direction)
        self.status = ''
        
        self.tstartpc = self.tlenpc = self.tendpc = self.qstartpc = self.qlenpc = self.qendpc = self.tscflen = self.qscflen = None
        if self.tscaffold in sequences:
            self.tscflen = len(sequences[self.tscaffold])
        if self.qscaffold in sequences:
            self.qscflen = len(sequences[self.qscaffold])

        if self.tscflen and self.qscflen:
            self.set_pc(1, self.tscflen, 1, self.qscflen)

        self.pbscaffold = None
        if args:
            self.pbscaffold = args[0]
            self.pbtstart, self.pbtlen, self.pbtend, self.pbtdir   = args[1], args[2], args[3], args[4]
            self.pbqstart, self.pbqlen, self.pbqend, self.pbqdir   = args[5], args[6], args[7], args[8]
            self.direction = self.pbtdir * self.pbqdir

    def set_pc(self, t_raftstart, t_raftend, q_raftstart, q_raftend):
        tlen = t_raftend - t_raftstart + 1
        qlen = q_raftend - q_raftstart + 1
        self.tstartpc = (self.tstart - t_raftstart + 1) / tlen * 100
        self.tlenpc   =  self.tlen                      / tlen * 100
        self.tendpc   = (self.tend   - t_raftstart)     / tlen * 100
        self.qstartpc = (self.qstart - q_raftstart + 1) / qlen * 100
        self.qlenpc   =  self.qlen   / qlen * 100
        self.qendpc   = (self.qend   - q_raftstart + 1) / qlen * 100
        
    def set_status(self, t_start=None, t_end=None, q_start=None, q_end=None):
        self.status = ''
        if self.tscflen and not (t_start and t_end):
            t_start = 1
            t_end = self.tscflen
        if self.qscflen and not (q_start and q_end):
            q_start = 1
            q_end = self.qscflen

        if t_start and q_start:
            t_start, t_end = order(t_start, t_end)
            q_start, q_end = order(q_start, q_end)
            
            self.set_pc(t_start, t_end, q_start, q_end)
            
            t_before = self.tstart - t_start
            t_after  = t_end - self.tend
            q_before = self.qstart - q_start
            q_after  = q_end - self.qend
            
            if self.pbscaffold is None:
                if self.direction == -1:
                    q_before, q_after = q_after, q_before
                
                if t_before < q_before and t_after < q_after:
                    self.status = 'within'
                    return

                if self.direction == 1:
                    if self.tstartpc < 10 and self.qendpc > 90:
                        self.status = 'connect_bf_af'
                    elif self.tendpc > 90 and self.qstartpc < 10:
                        self.status = 'connect_af_bf'
                elif self.direction == -1:
                    if self.tstartpc < 10 and self.qstartpc < 10:
                        self.status = 'connect_br_af'
                    elif self.tendpc > 90 and self.qendpc > 90:
                        self.status = 'connect_af_br'
                return
            else:
                if self.pbtstart < self.pbtend < self.pbqstart < self.pbqend:
                    if self.pbtdir == 1 and self.pbqdir == 1:
                        self.status = 'connect_af_bf'
                    elif self.pbtdir == -1 and self.pbqdir == 1:
                        self.status = 'connect_ar_bf'
                    elif self.pbtdir == 1 and self.pbqdir == -1:
                        self.status = 'connect_af_br'
                    else:
                        self.status = 'connect_ar_br'
                    return
                elif self.pbqstart < self.pbqend < self.pbtstart < self.pbtend:
                    if self.pbtdir == 1 and self.pbqdir == 1:
                        self.status = 'connect_bf_af'
                    elif self.pbtdir == -1 and self.pbqdir == 1:
                        self.status = 'connect_bf_ar'
                    elif self.pbtdir == 1 and self.pbqdir == -1:
                        self.status = 'connect_br_af'
                    else:
                        self.status = 'connect_br_ar'
                    return
                
        self.status = ''

    def __repr__(self):
        output = "{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{}\t{}\t{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{}\t{}".format(
                self.tscaffold, self.tstart, self.tlen, self.tend, self.tstartpc, self.tlenpc, self.tendpc,
                self.qscaffold, self.qstart, self.qlen, self.qend, self.qstartpc, self.qlenpc, self.qendpc, self.direction, self.status)
        if self.pbscaffold:
            output += "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.pbscaffold, self.pbtstart, self.pbtlen, self.pbtend, self.pbtdir, self.pbqstart, self.pbqlen, self.pbqend, self.pbqdir)
        return output

def order(start, end):
    if start < end:
        return start, end
    else:
        return end, start

if __name__ == '__main__':
    pass