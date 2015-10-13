#!/usr/bin/python3 -u

import os
import sys
import pathlib
import argparse
import glob
import re
import sqlite3 as sql
from pprint import pprint
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from collections import defaultdict
from operator import itemgetter

class Part():
    def __init__(self, oldname, oldstart, oldend, newname, newstart, newend, strand, parttype, comment="", chromosome="", cm=""):
        self.oldname = oldname
        self.oldstart = None if oldstart is None else int(oldstart)
        self.oldend = None if oldend is None else int(oldend)
        self.newname = newname
        self.newstart = int(newstart)
        self.newend = int(newend)
        self.strand = int(strand)
        self.parttype = parttype
        self.comment = comment
        self.chromosome = chromosome
        self.cm = cm

    @property
    def length(self):
        return self.oldend - self.oldstart + 1
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.oldname, self.oldstart, self.oldend, self.newname,
                                                               self.newstart, self.newend, self.strand, self.parttype,
                                                               self.chromosome, self.cm)

class Comment:
    def __init__(self, name, oldstart, oldend, start, end, slice_start, slice_end, strand):
        self.name = name
        self.oldstart = oldstart
        self.oldend = oldend
        self.start = start
        self.end = end
        self.slice_start = slice_start
        self.slice_end = slice_end
        self.strand = strand

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.name, self.oldstart, self.oldend, self.start, self.end, self.slice_start, self.slice_end, self.strand)

class MapPart:
    def __init__(self, scaffold, start, end, length, chromosome, cm, comment=""):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.length = length
        self.chromosome =chromosome
        self.cm = cm
        self.comment = comment
    
    def __repr__(self):
        return '{}:{}-{} {}bp {} {}'.format(self.scaffold, self.start, self.end, self.length, self.chromosome, self.cm)


class MapScaffold:
    def __init__(self):
        self.mapparts = []
    
    def append(self, scaffold, start, end, length, chromosome, cm):
        self.mapparts.append(MapPart(scaffold, start, end, length, chromosome, cm))


def load_transfers(mergedgenome):
    genome = defaultdict(list)
    try:
        with open(mergedgenome, 'r') as g:
            for line in g:
                oldname, oldstart, oldend, newname, newstart, newend, strand, typename, *args = line.rstrip().split('\t')
                if typename in ["full", "part"]:
                    typename = "retained"
                genome[oldname].append(Part(oldname, oldstart, oldend, newname, newstart, newend, strand, typename))
    except IOError:
        print("Failed to load genome file {}".format(mergedgenome))
        sys.exit()
    
    return genome

def load_linkage_map(database, errorarg, genome):
    conn_in, ci = open_input_database(database)
    
    errors = load_errors(errorarg)

    linkage_map = {}
    for chromosome, cm, scaffold, start, end, length, *args in ci.execute('select * from scaffold_map order by scaffold, start'):

        if scaffold in errors and start in errors[scaffold]:
            chromosome = 0
            cm = -1

        if not scaffold in linkage_map:
            linkage_map[scaffold] = MapScaffold()

        linkage_map[scaffold].append(scaffold, start, end, length, chromosome, cm)
    
    for scaffold in genome:
        if scaffold not in linkage_map:
            linkage_map[scaffold] = MapScaffold()
            for part in genome[scaffold]:
                if part.parttype in ['active', 'retained']:
                    linkage_map[scaffold].append(scaffold, part.oldstart, part.oldend, part.length, 0, -1)

    return linkage_map, ci

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


def load_errors(errorarg):
    errors = defaultdict(lambda:defaultdict(int))
    if not errorarg:
        return errors

    try:
        if not os.path.isfile(errorarg):
            raise IOError
        with open(errorarg, 'r') as e:
            for line in e:
                scaffold, start = line.rstrip().split('\t')
                start = int(start)
                errors[scaffold][start] = 0

        return errors

    except IOError:
        print("Can't open error file {}".format(errorarg))
        sys.exit(1)


def load_genome(fasta):
    try:
        with open(fasta, 'r') as f:
            sequences = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
            for scaffold in sequences:
                sequences[scaffold].seq = sequences[scaffold].seq.upper()

    except IOError:
        print("Can't load genome from file {}!".format(fasta))
        sys.exit()

    return sequences
    

class GFF:
    def __init__(self, *args):
        if len(args) == 1:
            self.scaffold, self.source, self.featuretype, self.start, self.end, self.score, self.strand, self.phase, self.attributes = args[0].rstrip().split('\t')
        else:
            self.scaffold    = args[0]
            self.source      = args[1]
            self.featuretype = args[2]
            self.start       = args[3]
            self.end         = args[4]
            self.score       = args[5]
            self.strand      = args[6]
            self.phase       = args[7]
            self.attributes  = args[8]
        self.start = int(self.start)
        self.end = int(self.end)
        self.genename = get_gene_name(self)
        if self.attributes.endswith(';'):
            self.attributes = self.attributes[:-1]
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.scaffold, self.source, self.featuretype, self.start, self.end, self.score, self.strand, self.phase, self.attributes)

def load_crossmap(crossmapout):
    genes = defaultdict(list)
    gene_exon_positions = defaultdict(lambda: defaultdict(tuple))
    try:
        with open(crossmapout) as c:
            for line in c:
                cm = None
                if 'fail' in line:
                    featureline = line.split('\tfail')[0]
                if '->' in line:
                    featureline, cmline = line.split('\t->\t')
                    cm = GFF(cmline)
                    cm.attributes += ';Note=CrossMap'
                feature=GFF(featureline)
                feature.crossmap = cm
                gene_id = get_gene_attribute(feature, "ID")
                if feature.featuretype == 'exon':
                    gene_exon_positions[feature.genename][gene_id] = (feature.start, feature.end)
                if feature.featuretype == 'CDS':
                    for exon in gene_exon_positions[feature.genename]:
                        e = gene_exon_positions[feature.genename][exon]
                        if e[0] <= feature.start <= e[1] and e[0] <= feature.end <= e[1]:
                            gene_id = exon + "_CDS"
                    
                    if feature.crossmap:    # Lose CrossMap features where the CDSs aren't the same length
                        cds_len = feature.end - feature.start + 1
                        cm_len  = feature.crossmap.end - feature.crossmap.start + 1
                        if cds_len != cm_len:
                            feature.crossmap = None
                if gene_id is None:
                    print("No gene id for CDS found", feature, end="")
                feature.id = gene_id

                genes[feature.genename].append(feature)
                
    except IOError:
        print("Failed to load CrossMap output {}".format(crossmapout))
        sys.exit()

    return genes
    
def load_gff(gff):
    genes = defaultdict(list)
    gene_exon_positions = defaultdict(lambda: defaultdict(tuple))
    try:
        with open(gff) as g:
            for line in g:
                if line.startswith('#') or 'contig' in line:
                    continue
                feature = GFF(line)
                gene_id = get_gene_attribute(feature, "ID")
                if feature.featuretype == 'exon':
                    gene_exon_positions[feature.genename][gene_id] = (feature.start, feature.end)
                if feature.featuretype == 'CDS':
                    for exon in gene_exon_positions[feature.genename]:
                        e = gene_exon_positions[feature.genename][exon]
                        if e[0] <= feature.start <= e[1] and e[0] <= feature.end <= e[1]:
                            gene_id = exon + "_CDS"
                if gene_id is None:
                    print("No gene id for CDS found", feature, end="")
                feature.id = gene_id
                genes[feature.genename].append(feature)
                
    except IOError:
        print("Failed to load GFF file {}".format(gff))
        sys.exit()
    
    return genes

def get_gene_name(feature):
    genename = None
    if feature.featuretype in ['gene', 'pseudogenic_region']:
        genename = re.search(r"ID=(.+?)(;|$)", feature.attributes).group(1)
    else:
        genename = re.search(r"Parent=(.+?)([;-]|$)", feature.attributes).group(1)
    return genename

def get_gene_attribute(feature, attribute="ID"):
    gene_ids = re.search(r'{}=(.+?)(;|$)'.format(attribute), feature.attributes)
    if gene_ids:
        gene_id = gene_ids.group(1)
    else:
        gene_id = None
    return gene_id


def open_output_database(output):
    try:
        conn = sql.connect(output)
        db = conn.cursor()
        db.execute('drop table if exists scaffold_map')
        db.execute('''create table scaffold_map
                     (chromosome integer, cm real, scaffold text, start integer, end integer, length integer, parttype text, comments text)''')

        cursor = conn.cursor()
        return conn, cursor
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)


def get_parts(scaffold, start, end, genome):

    parts = []
    haps = []
    for part in genome[scaffold]:
        if part.oldstart > end or part.oldend < start:
            continue
        slice_start = max(part.oldstart, start)
        slice_end = min(part.oldend, end)
        slice_length = slice_end - slice_start + 1
        
        start_offset = slice_start - part.oldstart
        end_offset = part.oldend - slice_end
        if part.newstart is None:
            newstart = newend = None
        elif part.strand == 1 or part.strand == 0:
            newstart = part.newstart + start_offset
            newend = part.newend - end_offset
        elif part.strand == -1:
            newstart = part.newstart + end_offset
            newend = part.newend - start_offset

        newpart = Part(part.newname, newstart, newend, part.oldname, part.oldstart, part.oldend, part.strand, part.parttype,
                       Comment(part.oldname, part.oldstart, part.oldend, start, end, slice_start, slice_end, part.strand))
        if part.parttype in ['haplotype', 'gap']:
            haps.append(newpart)
        else:
            parts.append(newpart)
    return parts, haps

def transfer_chromosome_map(ci, co):
    
    table = ci.execute('select * from sqlite_master where name="chromosome_map" and type="table"')
    if len(table.fetchall()) != 1:
        return
    try:
        co.execute('drop table if exists chromosome_map')
        co.execute('''create table chromosome_map
                     (chromosome integer, print text, cm real, original text, clean text, length integer)''')
        
        chrommap = ci.execute('select * from chromosome_map')
        for (chromosome, chromprint, cm, original, clean, length) in chrommap.fetchall():
            co.execute('insert into chromosome_map values (?,?,?,?,?,?)', [chromosome, chromprint, cm, original, clean, length])
    except sql.Error as sqle:
        print(sqle)
        sys.exit(1)

def write_new_map(linkage_map, genome, output, ci, collapse):
    conn_out, co = open_output_database(output + "_map.db")

    transfer_chromosome_map(ci, co)

    new_map = []

    genomelength = 0
    for scaffold in linkage_map:
        for part in linkage_map[scaffold].mapparts:
            new_parts, haps = get_parts(scaffold, part.start, part.end, genome)
            for np in new_parts:
                np.chromosome = part.chromosome
                np.cm = part.cm
                new_map.append(np)

    new_map.sort(key=lambda x: (x.oldname, x.oldstart))

    if collapse:
        collapse_map(new_map)

    for mp in new_map:
        co.execute('insert into scaffold_map values (?,?,?,?,?,?,?,?)', [mp.chromosome, mp.cm, mp.oldname, mp.oldstart, mp.oldend, mp.oldend-mp.oldstart+1, mp.parttype, repr(mp.comment)])

    conn_out.commit()
    conn_out.close()

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
                        blockmerge(mpa, mpb)
                    for mi in reversed(m):
                        del mapparts[mi]
                break
            j += 1
        
        if not merge:
            i += 1

    fill_maternal_only_cms(mapparts)
    recollapse(mapparts) # Collapses new filled maternal blocks and blocks with no linkage information

def recollapse(mapparts):
    i = 0
    
    while i < len(mapparts) - 1:
        mpi = mapparts[i]
        j = i + 1
        if j == len(mapparts):
            break
        while (mpi.chromosome == mapparts[j].chromosome and mpi.cm == mapparts[j].cm and
               mpi.oldname == mapparts[j].oldname and mpi.comment.name == mapparts[j].comment.name and
               mpi.oldend == mapparts[j].oldstart-1):
            blockmerge(mpi, mapparts[j])
            del mapparts[j]
        i += 1
    
def fill_maternal_only_cms(mapparts):
    i = 0
    for i in range(len(mapparts)-2):
        for o in ([i, i+1, i+2], [i+2, i+1, i]):
            mp1, mp2, mp3 = mapparts[o[0]], mapparts[o[1]], mapparts[o[2]]
            if mp1.oldname != mp2.oldname or mp1.oldname != mp3.oldname or mp2.oldname != mp3.oldname:
                continue
            if (mp1.chromosome == 0 and mp1.cm == -1 and
                mp2.chromosome != 0 and mp2.cm == -1 and
                mp3.chromosome != 0 and mp3.cm != -1):
                mp2.cm = mp3.cm


def get_merge_groups(mapparts, i, parts):
    merge_groups = [[i]]
    merge_group_i = 0
    current_oldname = mapparts[i].oldname
    current_commentname = mapparts[i].comment.name
    
    for p in parts:
        mapparts[p].chromosome = mapparts[i].chromosome
        mapparts[p].cm = mapparts[i].cm
        if mapparts[p].oldname != current_oldname or mapparts[p].comment.name != current_commentname:
            merge_group_i += 1
            merge_groups.append([])
            current_oldname = mapparts[p].oldname
            current_commentname = mapparts[p].comment.name
        merge_groups[merge_group_i].append(p)
    return merge_groups

def blockmerge(a,b):
    a.oldend = b.oldend
    if a.comment.strand == 1:
        a.comment.oldend = b.comment.oldend
        a.comment.end = b.comment.end
        a.comment.slice_end = b.comment.slice_end
    elif a.comment.strand == -1:
        a.comment.oldstart = b.comment.oldstart
        a.comment.start = b.comment.start
        a.comment.slice_start = b.comment.slice_start

def collapse_removed(genome):
    for scaffold in sorted(genome):
        rejecttypes = ['haplotype', 'removed']
        i = 0
        while i < len(genome[scaffold])-1:
            parts = sorted(genome[scaffold], key = lambda x: x.oldstart)
            parti = parts[i]
            partj = parts[i+1]
            if parti.oldname == partj.oldname and parti.oldend+1 == partj.oldstart:
                if parti.parttype in rejecttypes and partj.parttype in rejecttypes:
                    parti.oldend = partj.oldend
                    parti.newname = parti.newstart = parti.newend = None
                    parti.parttype = 'removed'
                    for k, part in enumerate(genome[scaffold]):
                        if part == partj:
                            del genome[scaffold][k]
                    continue
            i += 1

def make_output_feature(outstart, outend, feature, new_parts):
    if outstart is not None and outend is not None:
        if outstart > outend:
            outstart, outend = outend, outstart
        strand = feature.strand
        if new_parts[0].strand == -1:
            if feature.strand == '+':
                strand = '-'
            elif feature.strand == '-':
                strand = '+'
        return GFF(new_parts[0].oldname, feature.source, feature.featuretype, outstart, outend, feature.score, strand, feature.phase, feature.attributes)
    else:
        return None

def transfer_gff_feature(feature, genome, crossmap=False):
    output_feature = outstart = outend = status = None

    if crossmap and feature.crossmap:
        output_feature = GFF(feature.crossmap.scaffold, feature.crossmap.source, feature.crossmap.featuretype, feature.crossmap.start, feature.crossmap.end, feature.crossmap.score, feature.crossmap.strand, feature.crossmap.phase, feature.crossmap.attributes)
        status = 'crossmap'
        return output_feature, status

    new_parts, haps = get_parts(feature.scaffold, feature.start, feature.end, genome)

    if not new_parts:
        if not haps:
            status = 'missing'
        else:
            status = haps[0].parttype
    elif len(new_parts) == 1 and not haps:
        np = new_parts[0]
        if 'removed' not in np.parttype:
            outstart, outend = new_parts[0].oldstart, new_parts[0].oldend
            status = 'ok'
        else:
            status = 'removed'
    else:
        scaffolds = {}
        for part in new_parts:
            if 'removed' in part.parttype:
                continue
            if part.newstart <= feature.start <= part.newend:
                outstart = part.oldstart if part.strand == 1 else part.oldend
                scaffolds[part.oldname] = 1
            if part.newstart <= feature.end <= part.newend:
                outend = part.oldend if part.strand == 1 else part.oldstart
                scaffolds[part.oldname] = 1
        if outstart is not None and outend is not None:
            status = 'ok'
        if len(scaffolds) > 1:
            outstart = outend = None
            status =  'multiscaffold'
        if feature.featuretype == 'CDS':    # Do not allow CDSs to span multiple parts, but try to find a CrossMap hit
            outstart = outend = None
            status = 'broken'
        
    output_feature = make_output_feature(outstart, outend, feature, new_parts)
    
    if output_feature is None and feature.crossmap:
        output_feature = GFF(feature.crossmap.scaffold, feature.crossmap.source, feature.crossmap.featuretype, feature.crossmap.start, feature.crossmap.end, feature.crossmap.score, feature.crossmap.strand, feature.crossmap.phase, feature.crossmap.attributes)
        status = 'crossmap'

    if output_feature is None:
        if status is None:
            status = 'broken'
        return GFF(feature.scaffold, feature.source, feature.featuretype, feature.start, feature.end, feature.score, feature.strand, feature.phase, feature.attributes), status
    else:
        return output_feature, status

def validate_gene(genefeatures):
    gene = None
    RNAs = []
    exons = defaultdict(list)
    cdses = defaultdict(list)
    
    for feature in genefeatures:
        if feature.featuretype == 'gene':
            gene = feature
        if 'RNA' in feature.featuretype:
            RNAs.append(feature)
            RNA_name = get_gene_attribute(feature, "ID")
        if feature.featuretype == 'exon':
            exons[RNA_name].append(feature)
        if feature.featuretype == 'CDS':
            cdses[RNA_name].append(feature)
            cds_has_exon = 0
            cds_parent = get_gene_attribute(feature, "Parent")
            for exon in exons[RNA_name]:
                if cds_parent == RNA_name and exon.start <= feature.start <= exon.end and exon.start <= feature.end <= exon.end:
                    cds_has_exon += 1
            if cds_has_exon == 0:
                return 'cds_not_in_exon'

    for RNA in cdses:
        mcds = cdses[RNA]
        orient = mcds[0].strand
        for i in range(len(mcds)-1):
            if orient == '+' and mcds[i].start > mcds[i+1].start or orient == '-' and mcds[i].start < mcds[i+1].start:
                return 'disordered_cds'

        for cds in cdses[RNA]:
            for cds2 in cdses[RNA]:
                if cds is cds2:
                    continue
                if cds.start <= cds2.start <= cds.end or cds.start <= cds2.end <= cds.end:
                    return 'cds_overlap'

    gene_start = 1000000000000000
    gene_end   = -1
    exon_strands = defaultdict(int)
    for exon_parent in exons:
        scaffold = exons[exon_parent][0].scaffold
        RNA_start = min(exon.start for exon in exons[exon_parent])
        RNA_end   = max(exon.end   for exon in exons[exon_parent])
        gene_start = min(RNA_start, gene_start)
        gene_end   = max(RNA_end, gene_end)
        for exon in exons[exon_parent]:
            exon_strands[exon.strand] += 1
    
        for RNA in RNAs:
            RNA_id = get_gene_attribute(RNA, "ID")
            if RNA_id != exon_parent:
                continue
            if RNA.scaffold != scaffold:
                RNA.scaffold = scaffold
            if RNA_start != RNA.start:
                RNA.start = RNA_start
            if RNA_end != RNA.end:
                RNA.end = RNA_end
    if gene:
        if gene.scaffold != scaffold:
            gene.scaffold = scaffold
        if gene_start != gene.start:
            gene.start = gene_start
        if gene_end != gene.end:
            gene.end = gene_end
        if len(exon_strands) == 1:
            strand = list(exon_strands)[0]
            gene.strand = strand
            for RNA in RNAs:
                RNA.strand = strand
        else:
            print("Too many strands for exons!")
            print(genefeatures)

    return 'ok'

class ValidRNA:
    def __init__(self, name, cdses, cds_lengths, start_ok, stop_ok, extra_stops, sequence):
        self.name = name
        self.cdses = cdses
        self.cds_lengths = cds_lengths
        self.start_ok = start_ok
        self.stop_ok = stop_ok
        self.extra_stops = extra_stops
        self.sequence = sequence
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.name, self.cdses, self.cds_lengths, self.start_ok, self.stop_ok, self.extra_stops)

def validate_protein(features, genome):
    stop_codons = ['TAA','TAG','TGA']
    
    valid_mRNAs = []
    mrnas = defaultdict(list)
    for feature in features:
        if feature.featuretype == 'mRNA':
            mrna_name = get_gene_attribute(feature, "ID")
        if feature.featuretype == 'CDS':
            mrnas[mrna_name].append(feature)

    for mrna in sorted(mrnas):
        seqrecord = None
        cds_lengths = []
        for cds in mrnas[mrna]:
            orientation = cds.strand
            cds_lengths.append(cds.end-cds.start+1)
            cds_seq = genome[cds.scaffold][(cds.start-1):cds.end]
            if cds.strand == '-':
                cds_seq = cds_seq.reverse_complement(id=True, name=True, description=True)
            if not seqrecord:
                seqrecord = cds_seq
            else:
                seqrecord += cds_seq
        
        start_ok = stop_ok = False
        extra_stops = 0
        if seqrecord[:3].seq == 'ATG':
            start_ok = True
        if seqrecord[-3:].seq in stop_codons:
            stop_ok = True
        extra_stops = 0
        for i in range(0, len(seqrecord)-3, 3):
            if seqrecord[i:i+3].seq in stop_codons:
                extra_stops += 1
        
        valid_mRNAs.append(ValidRNA(mrna, len(mrnas[mrna]),cds_lengths, start_ok, stop_ok, extra_stops, seqrecord.seq))
        
    return valid_mRNAs

def compare_proteins(genename, genes, output_features, oldseqs, newseqs):
    if oldseqs is None or newseqs is None:
        return 'ok'
    old_RNAs = validate_protein(genes[genename], oldseqs)
    new_RNAs = validate_protein(output_features, newseqs)
    if len(old_RNAs) != len(new_RNAs):
        return 'different_RNA_number'
    for i in range(len(old_RNAs)):
        o = old_RNAs[i]
        n = new_RNAs[i]
        if o.cdses != n.cdses:
            return 'different_CDS_number'
        for l in range(len(o.cds_lengths)):
            if o.cds_lengths[l] != n.cds_lengths[l]:
                return 'CDS_length_mismatch'
        if o.start_ok != n.start_ok:
            return 'bad_start_codon'
        if o.stop_ok != n.stop_ok:
            return 'bad_stop_codon'
        if o.extra_stops < n.extra_stops:
            return 'extra_stop_codons'
    return 'ok'

def get_output_features(gene, genetypes, genome, crossmap=False):
    output_features = []
    status = defaultdict(int)
    for feature in gene:
        output_feature, outstatus = transfer_gff_feature(feature, genome, crossmap)
        output_features.append(output_feature)
        if isinstance(output_feature, GFF) and output_feature.featuretype not in ['gene', 'mRNA']:
            status[outstatus] += 1

    statuses = sorted(status.keys())
    genetype = '_'.join(statuses)
    genetypes[genetype] += 1
    return output_features, genetype

def get_gene_status(genetype, genestats=None):
    if genetype == 'missing':
        gene_status = 'missing'
    elif 'crossmap-removed' in genetype:
        gene_status = 'removed'
    elif genetype in ['disordered_cds', 'cds_overlap', 'extra_stop_codons', 'bad_start_codon', 'bad_stop_codon', 'cds_not_in_exon']:
        gene_status = 'broken'
    elif 'broken' in genetype or 'multiscaffold' in genetype or 'missing' in genetype:
        gene_status = 'broken'
    elif 'ok' in genetype and ('haplotype' in genetype or 'gap' in genetype or 'removed' in genetype):
        gene_status = 'hap_broken'
    elif 'haplotype' in genetype or 'gap' in genetype or 'removed' in genetype:
        gene_status = 'removed'
    else:
        gene_status = 'ok'

    if genestats is not None:
        genestats[gene_status] += 1

    return gene_status

def try_validation(output_features, genename, genes, oldseqs, newseqs):
        
    scaffolds = defaultdict(int)
    for of in output_features:
        if isinstance(of, GFF) and of.featuretype not in ['gene', 'mRNA']:
            scaffolds[of.scaffold] = 1

    if len(scaffolds) > 1:
        return 'multiscaffold'

    validgene = validate_gene(output_features)
    if validgene != 'ok':
        return validgene

    validprotein = compare_proteins(genename, genes, output_features, oldseqs, newseqs)
    return validprotein


def transfer_gene(genename, genes, genetypes, genestats, genes_by_scaffold, gene_starts, genome, oldseqs, newseqs, brokenfile, removefile):
    
    output_features, genetype = get_output_features(genes[genename], genetypes, genome)
    genestatus = get_gene_status(genetype, genestats)

    if genestatus is not 'ok':
        write_rejected_gene(genename, genestatus, genes, brokenfile, removefile)
        return

    validstatus = try_validation(output_features, genename, genes, oldseqs, newseqs)

    if validstatus != 'ok':
        genetypes[genetype] -= 1
        output_features, genetype = get_output_features(genes[genename], genetypes, genome, True)
        validstatus = try_validation(output_features, genename, genes, oldseqs, newseqs)
        if validstatus != 'ok':
            reject_gene(validstatus, genename, genetype, genestatus, genetypes, genestats, genes, brokenfile, removefile)
            return
    
    for of in output_features:
        if 'gene' in of.featuretype or 'pseudogenic_region' in of.featuretype:
            gene_starts[genename] = of.start
        genes_by_scaffold[of.scaffold][genename].append(of)

    return


def reject_gene(reason, genename, genetype, genestatus, genetypes, genestats, genes, brokenfile, removefile):
    
    genestats[genestatus] -= 1
    
    newstatus = 'broken'
    if genetype == 'crossmap':
        newstatus = 'removed'
        reason = reason+'-crossmap-removed'

    genestats[newstatus] += 1
    genetypes[genetype] -= 1
    genetypes[reason] += 1

    write_rejected_gene(genename, newstatus, genes, brokenfile, removefile)

def write_rejected_gene(genename, genestatus, genes, brokenfile, removefile):
    for feature in genes[genename]:
        if genestatus in ['broken', 'hap_broken']:
            brokenfile.write(repr(feature))
        elif genestatus in ['missing', 'removed']:
            removefile.write(repr(feature))


def write_ok_gff(output, genes_by_scaffold, gene_starts, newseqs):
    gfffile = open(output + ".gff", 'w')
    gfffile.write("##gff-version 3\n")

    if newseqs:
        for scaffold in newseqs:
            gfffile.write("##sequence-region {} 1 {}\n".format(scaffold, len(newseqs[scaffold].seq)))
        
    for scaffold in sorted(genes_by_scaffold):
        for genename in sorted(genes_by_scaffold[scaffold], key = lambda x: gene_starts[x]):
            for feature in genes_by_scaffold[scaffold][genename]:
                gfffile.write(repr(feature))

    gfffile.close()

def write_new_gff(genes, genome, oldseqs, newseqs, output):

    genestats = defaultdict(int)
    genetypes = defaultdict(int)
    genes_by_scaffold = defaultdict(lambda: defaultdict(list))
    gene_starts = defaultdict(int)

    removefile = open(output + '_removed.gff', 'w')
    removefile.write("##gff-version 3\n")
    
    brokenfile = open(output + '_broken.gff', 'w')
    brokenfile.write("##gff-version 3\n")

    for genename in genes:
        transfer_gene(genename, genes, genetypes, genestats, genes_by_scaffold, gene_starts, genome, oldseqs, newseqs, brokenfile, removefile)

    removefile.close()
    brokenfile.close()

    write_ok_gff(output, genes_by_scaffold, gene_starts, newseqs)

    print_gene_stats(genestats, genetypes)
    

def print_gene_stats(genestats, genetypes):
    total = sum(genestats.values())
    print_stat("OK", genestats['ok'], total)
    print_stat("Removed", genestats['removed'], total)
    print_stat("Brokehap", genestats['hap_broken'], total)
    print_stat("Broken", genestats['broken'], total)
    print_stat("Missing", genestats['missing'], total)
    print("Total   genes:\t{}".format(total))

    type_genestats = defaultdict(int)
    for genetype in sorted(genetypes, key=genetypes.get, reverse=True):
        print("{}\t{}".format(genetypes[genetype],genetype))
        type_genestats[get_gene_status(genetype)] += genetypes[genetype]
    for stat in sorted(type_genestats, key=type_genestats.get, reverse=True):
        print("{}\t{}".format(type_genestats[stat], stat))

def print_stat(text, stat, total):
    if total > 0:
        summary = stat/total*100
    else:
        summary = 0
    print ('{:<8} genes:\t{:>5}\t{:5.2f} %'.format(text, stat, summary))


def get_args():
    parser = argparse.ArgumentParser(description='''Output new database containing old linkage information on new HaploMerger output
    
        -m mergedgenome
        -g gff
        -d database
        -e errors
        -o output
        -c crossmap
        -f fasta
        -n newgenomefasta
        -x collapseoff
        ''')

    parser.add_argument('-m', '--mergedgenome', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=False)
    parser.add_argument('-d', '--database', type=str, required=False)
    parser.add_argument('-e', '--errors', type=str, required=False)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-c', '--crossmap', type=str, required=False)
    parser.add_argument('-f', '--fasta', type=str, required=False)
    parser.add_argument('-n', '--newgenomefasta', type=str, required=False)
    parser.add_argument('-x', '--collapse', action='store_false', default=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    genome = load_transfers(args.mergedgenome)
    
    if args.database:
        linkage_map, ci = load_linkage_map(args.database, args.errors, genome)
        write_new_map(linkage_map, genome, args.output, ci, args.collapse)

    collapse_removed(genome)

    genes = None
    if args.gff:
        genes = load_gff(args.gff)

    crossmap = defaultdict(list)
    if args.crossmap:
        genes = load_crossmap(args.crossmap)

    if genes:
        oldseqs = newseqs = None
        if args.fasta:
            oldseqs = load_genome(args.fasta)
        if args.newgenomefasta:
            newseqs = load_genome(args.newgenomefasta)

        write_new_gff(genes, genome, oldseqs, newseqs, args.output)