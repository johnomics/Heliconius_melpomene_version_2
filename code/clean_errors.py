#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import glob
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict, namedtuple
from operator import itemgetter

class MapPart:
    def __init__(self, scaffold, start, end, chromosome, cm):
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.chromosome = int(chromosome)
        self.cm = float(cm)

    @property
    def length(self):
        return self.end - self.start + 1

    def __repr__(self):
        return '{}:{:>8}-{:>8} {:>8}bp\t{}\t{}'.format(self.scaffold, self.start, self.end, self.length, self.chromosome, self.cm)

def collapse_map(mapparts):
    i = 0
    while i < len(mapparts) - 1:
        mpi = mapparts[i]
        if mpi.chromosome == 0 or mpi.cm == -1:
            i += 1
            continue
        
        parts = []
        j = i + 1
        merge = False

        while j < len(mapparts):
            mpj = mapparts[j]
            parts.append(j)

            if mpi.chromosome == mpj.chromosome and mpj.cm == -1:
                j += 1
                continue

            # If j chromosome is real and chromosomes or cms don't match, skip
            if mpj.chromosome != 0 and mpj.cm != -1 and  (mpj.chromosome != mpi.chromosome or mpj.cm != mpi.cm):
                if mpi.chromosome == mpj.chromosome: # Even if cms don't match, if chromosomes match, fill them
                    for p in parts:
                        mapparts[p].chromosome = mpi.chromosome
                break

            # If chromosomes and cms match, fill intermediate blocks
            if mpj.chromosome == mpi.chromosome and mpj.cm == mpi.cm:
                for p in parts:
                    mapparts[p].chromosome = mpi.chromosome
                    mapparts[p].cm = mpi.cm
                
                # Split parts into groups belonging to the same scaffolds
                merge_groups = get_merge_groups(mapparts, i, parts)

                # Merge groups belonging to the same scaffolds
                for m in reversed(merge_groups):
                    if len(m) == 1:
                        continue
                    merge=True
                    mpa = mapparts[m[0]]
                    m.pop(0)
                    for mi in m:
                        mpb = mapparts[mi]
                        mpa.end = mpb.end
                    for mi in reversed(m):
                        del mapparts[mi]
                break
            j += 1
        
        if not merge:
            i += 1

    fill_maternal_only_cms(mapparts)
    recollapse(mapparts) # Collapses new filled maternal blocks and blocks with no linkage information

def get_merge_groups(mapparts, i, parts):
    merge_groups = [[i]]
    merge_group_i = 0
    current_scaffold = mapparts[i].scaffold
    
    for p in parts:
        mapparts[p].chromosome = mapparts[i].chromosome
        mapparts[p].cm = mapparts[i].cm
        if mapparts[p].scaffold != current_scaffold:
            merge_group_i += 1
            merge_groups.append([])
            current_scaffold = mapparts[p].scaffold
        merge_groups[merge_group_i].append(p)
    return merge_groups

def recollapse(mapparts):
    i = 0
    while i < len(mapparts) - 1:
        mpi = mapparts[i]
        j = i + 1
        while j < len(mapparts) and mpi.chromosome == mapparts[j].chromosome and mpi.cm == mapparts[j].cm:
            mpi.end = mapparts[j].end
            del mapparts[j]
        i += 1
    
def fill_maternal_only_cms(mapparts):
    i = 0
    for i in range(len(mapparts)-2):
        for o in ([i, i+1, i+2], [i+2, i+1, i]):
            mp1, mp2, mp3 = mapparts[o[0]], mapparts[o[1]], mapparts[o[2]]
            if (mp1.chromosome == 0 and mp1.cm == -1 and
                mp2.chromosome != 0 and mp2.cm == -1 and
                mp3.chromosome != 0 and mp3.cm != -1):
                mp2.cm = mp3.cm

class Stats:
    def __init__(self):
        self.blocks_loaded = 0
        self.bases_loaded = 0
        self.blocks_covered = 0
        self.bases_covered = 0
        self.error_blocks = 0
        self.error_bases = 0

def update_stats(s, length, chromosome):
    s.blocks_loaded += 1
    s.bases_loaded += length
    
    if chromosome != 0:
        s.blocks_covered += 1
        s.bases_covered += length
    
def clean_linkage_map(database, errorarg, reversearg, genome):

    conn_in, ci = open_input_database(database)
    
    errors = load_errors(errorarg)
    
    cms = load_chromosome_map(ci)
    s = Stats()

    genome = {}
    for chromosome, cm, scaffold, start, end, length, *args in ci.execute('select * from scaffold_map order by scaffold, start'):
        
        if scaffold in errors and start in errors[scaffold]:
            if errors[scaffold][start] == "D":
                del errors[scaffold][start]
                continue
            elif errors[scaffold][start] == "N":
                chromosome = 0
                cm = -1
                s.error_blocks += 1
                s.error_bases += length
                del errors[scaffold][start]

        update_stats(s, length, chromosome)

        if scaffold not in genome:
            genome[scaffold] = []
        if start not in errors[scaffold]:
            genome[scaffold].append(MapPart(scaffold, start, end, chromosome, cm))
    
    for scaffold in sorted(errors):
        for start in sorted(errors[scaffold]):
            end, chromosome, cm = errors[scaffold][start].split(',')
            genome[scaffold].append(MapPart(scaffold, start, end, chromosome, cm))
            update_stats(s, int(end)-int(start)+1, chromosome)
            
            chromosome, cm = int(chromosome), float(cm)
            if cm not in cms[chromosome]:
                cms[chromosome][cm] = cm

    cms = do_reverse(cms, reversearg)

    for scaffold in genome:
        genome[scaffold].sort(key = lambda x:x.start)
        collapse_map(genome[scaffold])
        
        for p in genome[scaffold]:
            comment = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(p.scaffold, p.start, p.end, p.start, p.end, p.start, p.end, 1)
            co.execute('insert into scaffold_map values (?,?,?,?,?,?,?,?)', 
                        [p.chromosome, cms[p.chromosome][p.cm], p.scaffold, p.start, p.end, p.length, "active", comment])

    print("Loaded  {} blocks and {} bases".format(s.blocks_loaded, s.bases_loaded))
    print("Covered {} blocks and {} bases".format(s.blocks_covered, s.bases_covered))
    print("Errors  {} blocks and {} bases".format(s.error_blocks, s.error_bases))
    return conn_in, ci, cms

def open_input_database(database):
    try:
        if not os.path.exists(database):
            raise IOError
        conn = sql.connect(database)
        cursor = conn.cursor()
        return conn, cursor
    except IOError:
        print("Can't open database {}".format(database))
        sys.exit(1)
    except sql.Error:
        print("SQL error")
        sys.exit(1)

def load_reverse_file(reversearg):

    reverses = {}
    if not reversearg:
        return reverses

    try:
        if not os.path.isfile(reversearg):
            raise IOError
        with open(reversearg, 'r') as r:
            for line in r:
                chromosome = line.rstrip()
                reverses[int(chromosome)]=True
    except IOError:
        print("Can't open error file {}".format(reversearg))
        sys.exit(1)

    return reverses
    
def do_reverse(cms, reversearg):

    reverses = load_reverse_file(reversearg)
    
    for chrom in cms:
        if chrom in reverses:
            maxcm = max(cms[chrom])
            for cm in cms[chrom]:
                if cm == -1.0:
                    continue
                revcm = format(maxcm - cm, '.3f')
                cms[chrom][cm] = revcm

    return cms


def load_chromosome_map(ci):

    query = ci.execute('select distinct chromosome from chromosome_map')
    chromosomes = [x[0] for x in query.fetchall()]
    
    cms = defaultdict(lambda: defaultdict(int))
    cms[0][-1.0] = -1.0
    for chrom in chromosomes:
        query = ci.execute('select distinct cm from chromosome_map where chromosome == {}'.format(chrom))
        cms[chrom][-1.0] = -1.0
        for cm in [x[0] for x in query.fetchall()]:
            cms[chrom][cm] = cm

    return cms


def load_errors(errorarg):
    errors = defaultdict(lambda:defaultdict(int))
    if not errorarg:
        return errors

    try:
        if not os.path.isfile(errorarg):
            raise IOError
        with open(errorarg, 'r') as e:
            for line in e:
                scaffold, start, action = line.rstrip().split('\t')
                start = int(start)
                errors[scaffold][start] = action

        return errors

    except IOError:
        print("Can't open error file {}".format(errorarg))
        sys.exit(1)

def open_output_database(output):
    try:
        conn = sql.connect(output)
        db = conn.cursor()
        db.execute('drop table if exists scaffold_map')
        db.execute('''create table scaffold_map
                     (chromosome integer, cm real, scaffold text, start integer, end integer, length integer, parttype text, comment text)''')

        cursor = conn.cursor()
        return conn, cursor
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)

def transfer_chromosome_map(ci, co, cms):
    try:
        co.execute('drop table if exists chromosome_map')
        co.execute('''create table chromosome_map
                     (chromosome integer, print text, cm real, original text, clean text, length integer)''')
        
        chrommap = ci.execute('select * from chromosome_map')
        for (chromosome, chromprint, cm, original, clean, length) in chrommap.fetchall():
            cm = cms[chromosome][cm] # Reverse chromosomes where necessary
            co.execute('insert into chromosome_map values (?,?,?,?,?,?)', [chromosome, chromprint, cm, original, clean, length])
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -d database
        -e errors
        -o output
        -r reverses
        ''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-r', '--reverses', type=str, required=False)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    conn_out, co = open_output_database(args.output)

    conn_in, ci, cms = clean_linkage_map(args.database, args.errors, args.reverses, co)

    transfer_chromosome_map(ci, co, cms)

    conn_in.close()
    conn_out.commit()
    conn_out.close()
