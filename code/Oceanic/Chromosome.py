#!/usr/bin/env python3

from . import GenomeData as gd
from . import Pool as pl
from . import Raft as r
from . import Stats as s
from . import Mergers as merge

import multiprocessing
import sys
from os.path import isfile
from math import log
import sqlite3 as sql
from collections import defaultdict

class Chromosome:
    def __init__(self, name, genome):
        self.name = str(name)
        self.markers = {}
        self.new_cms, self.old_cms = self.calculate_real_cms(self.name, genome.db)
        self.set_markers(genome.db)
        self.pools = []
        self.scaffold_map = []
        self.chromosome_map = []
        self.agp = []
        self.unmapped_start = 1
        self.unmapped_part  = 1
        self.genome = genome
        self.mapped_blocks, self.mapped_blocks_length, self.placed_blocks, self.placed_blocks_length = self.set_blocks()
        self.scaffold_count = 1
        
    def __repr__(self):
        return('\n'.join([self.name, repr(self.pools)]))

    def __iter__(self):
        return iter(self.pools)

    def threadstart(self, args):
        conn = self.open_database(args.database)
        self.db = conn.cursor()

    def open_database(self, dbfile):
        try:
            if isfile(dbfile):
                conn = sql.connect(dbfile)
            else:
                raise IOError
        except:
            print("Can't open database {}!".format(dbfile))
            sys.exit()
    
        return conn


    @property
    def stats(self):
        stats = s.Stats(self.name)
        stats.pool_num = len(self.pools)
        for pool in self.pools:
            gt = pool.pooltype
            for raft in pool:
                stats.length += raft.length
                stats.log_num += len(raft.logs)
                stats.pool_stats[gt].length += raft.length
                stats.log_count[gt] += len(raft.logs)
                stats.raft_count[gt] += 1
            stats.pool_stats[gt].count += 1
        
        return stats

    def set_markers(self, db):

        prev_cm = -1
        next_cm = -1
        for cm, in db.execute("select distinct cm from scaffold_map where chromosome={} order by cm".format(self.name)):

            if cm == -1:
                self.add_marker(-1, -1, -1)
                continue
            
            newcm = self.new_cms[cm]
            self.add_marker(newcm, prev_cm=prev_cm)
            if prev_cm != -1:
                self.update_marker(prev_cm, next_cm = newcm)
        
            prev_cm = newcm


    def calculate_real_cms(self, chromosome, db):
        cm_marker = {}
        for cm, pattern in db.execute("select cm, clean from chromosome_map where chromosome={} group by cm order by cm".format(chromosome)):
            marker = ''.join([x for x in pattern if x in ['A','B','H']])
            cm_marker[cm] = marker
        
        new_cms = {}
        new_cms[0.0] = 0.0
        new_cms[-1]  = -1
        old_cms = {}
        old_cms[0.0] = 0.0
        old_cms[-1]  = -1
        cmlist = list(sorted(cm_marker))
        cumul_cm = 0
        for i in range(1,len(cmlist)):
            a = cm_marker[cmlist[i-1]]
            b = cm_marker[cmlist[i]]
            diff = sum (a[j] != b[j] for j in range(len(a)))
            ind = len(a)
            newcm = 25 * log((1+2*diff/ind)/(1-2*diff/ind)) # Kosambi mapping function
            cumul_cm += newcm
            outcm = float('{:.3f}'.format(cumul_cm))
            new_cms[cmlist[i]] = outcm
            old_cms[outcm] = cmlist[i]
        return new_cms, old_cms

    def set_blocks(self):
        mapped_blocks = mapped_blocks_length = placed_blocks = placed_blocks_length = 0

        scaffold_cms = defaultdict(lambda: defaultdict(int))
        for newcm in sorted(self.markers):
            cm_blocks = pl.Pool(self)
            cm_block_id = 1
            statement = "select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(self.name, self.old_cms[newcm])
            for scaffold, start, end, length in self.genome.db.execute(statement):
                if scaffold not in self.genome.sequences:
                    continue
                if (scaffold, start) in self.genome.errors:
                    continue
                mapped_blocks += 1
                mapped_blocks_length += length
                self.genome.blocks[scaffold][start].add_marker(self.name,newcm)

                scaffold_cms[scaffold][newcm] = 1

                if newcm != -1:
                    placed_blocks += 1
                    placed_blocks_length += length
                
                cm_blocks.add(r.Raft(cm_block_id, scaffold, start, self))
                cm_block_id += 1

            if newcm != -1:
                self.pools.append(cm_blocks)

        return mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length


    def add_marker(self, cm, prev_cm=-1, next_cm=-1):
        self.markers[cm] = Marker(cm, prev_cm, next_cm)
        
    def update_marker(self, cm, prev_cm = -1, next_cm = -1):
        if cm not in self.markers:
            self.add_marker(cm, prev_cm, next_cm)

        self.markers[cm].update_previous(prev_cm)
        self.markers[cm].update_next(next_cm)

    def write(self):
        scaffolds = []
        
        if self.genome.revised:
            self.revised_conn = self.open_database(self.genome.revised+".db")
            self.revised_db = self.revised_conn.cursor()
        
        for pool in self:
            scaffolds += pool.write()
        
        self.revised_conn.commit()
        
        self.make_summaries()
        
        return scaffolds

    def update_tsv(self):
        parts = {}
        for pool in self.pools:
            for raft in pool.rafts:
                edges = sorted([raft.start, raft.end])
                for edge in edges:
                    i = 0
                    while i < len(self.genome.newparts[raft.scaffold]):
                        part = self.genome.newparts[raft.scaffold][i]
                        oldname = part.oldname # the Part referenced by part changes later, so copy name here
                        if part.parttype not in  ['haplotype'] and part.newstart < edge < part.newend:
                            if edge == edges[0]:
                                edge = edge - 1
                            offset = edge - part.newstart
                            if part.strand == '1':
                                oldedge = part.oldstart + offset
                                part1 = gd.OrigPart("\t".join([oldname, str(part.oldstart), str(oldedge), part.newname, str(part.newstart), str(edge), part.strand, part.parttype]))
                                part2 = gd.OrigPart("\t".join([oldname, str(oldedge+1), str(part.oldend), part.newname, str(edge+1), str(part.newend), part.strand, part.parttype]))
                            else:
                                oldedge = part.oldend - offset - 1
                                part1 = gd.OrigPart("\t".join([oldname, str(oldedge+1), str(part.oldend), part.newname, str(part.newstart), str(edge), part.strand, part.parttype]))
                                part2 = gd.OrigPart("\t".join([oldname, str(part.oldstart), str(oldedge), part.newname, str(edge+1), str(part.newend), part.strand, part.parttype]))
                            origi = -1
                            for j, oldpart in enumerate(self.genome.origparts[oldname]):
                                if oldpart.oldname == oldname and oldpart.oldstart == part.oldstart:
                                    origi = j
                                    break
                            self.genome.newparts[raft.scaffold][i] = part1
                            self.genome.newparts[raft.scaffold].insert(i+1, part2)
                            if origi != -1:
                                self.genome.origparts[oldname].insert(origi+1, part2)
                                self.genome.origparts[oldname][origi] = part1
                            else:
                                print("Can't find old part in original TSV: {}\n{}\n{}".format(raft.scaffold, part1, part2))
                            break
                        i += 1

    def discard_within(self):
        within_markers = []
        for pool in self.pools:
            for raft in pool:
                markers = raft.marker_chain
                if markers:
                    within = markers[1:(len(markers)-1)]
                    within_markers += within
        within_markers = set(within_markers)
        
        for pool in self.pools:
            for raft in pool:
                markers = raft.marker_chain
                if not markers:
                    continue
                remove = sum((1 for marker in markers if marker in within_markers))
                if remove == len(markers):
                    for pool_match in self.pools:
                        for raft_match in pool_match:
                            if raft_match is raft:
                                continue
                            match_markers = raft_match.marker_chain
                            if match_markers and markers[0] in match_markers and raft.scaffold != raft_match.scaffold:
                                raft.discard("within_removed")
            pool.cleanup()
        self.remove_empty_pools()

    def remove_mapped_offcuts(self):
        for i, pool in enumerate(self.pools):
            for raft in pool:
                remove = False
                if raft.offcuts:
                    for offcut in raft.offcuts:
                        if offcut in pool.scaffolds or i > 0 and offcut in self.pools[i-1].scaffolds or i < len(self.pools)-1 and offcut in self.pools[i+1].scaffolds:
                            remove = True
                if remove:
                    raft.discard("offcut_mapped_removed")
            pool.cleanup()
        self.remove_empty_pools()

    def remove_contained_nodes(self):
        for i, pool in enumerate(self.pools):
            for raft in pool:
                remove = False
                for scaffold in raft.scaffolds:
                    if scaffold in self.genome.nodes:
                        for qscaffold in self.genome.nodes[scaffold]:
                            if qscaffold in pool.scaffolds or (i > 0 and qscaffold in self.pools[i-1].scaffolds) or (i < len(self.pools)-1 and qscaffold in self.pools[i+1].scaffolds):
                                for node in self.genome.nodes[scaffold][qscaffold]:
                                    node.set_status(raft.start, raft.end)
                                    if node.status == 'within' and not self.genome.sequences[scaffold].features and (len(raft) >= 10000 and node.tlenpc > 50 or len(raft) < 10000 and node.tlenpc > 25):
                                        remove = True
                if remove:
                    raft.discard("within_node_removed")
            pool.cleanup()
        self.remove_empty_pools()

    def run_merger(self, mergeclass):
        if type(mergeclass) is not list:
            mergeclass = [mergeclass]
        for pool in self.pools:
            for mc in mergeclass:
                pool.assemble(pool, mc)

        self.connect(mergeclass)

    def assemble(self, args):

        self.threadstart(args)

        self.run_merger(merge.MarkerMerge)

        for pool in self.pools:
            pool.extend()

        self.update_tsv()

        self.discard_within()

        self.remove_mapped_offcuts()

        self.remove_contained_nodes()
        
        self.run_merger(merge.NodeMerge)
        self.remove_empty_pools()

        self.connect(merge.OKMerge)
        self.remove_empty_pools()

        print(self)
        print(self.stats)


    def make_summaries(self):
        # Data frame for graphical output
        chrom_sbs = []
        chrstart = 1
        chrpart = 1
        for i, pool in enumerate(self.pools):
            k = 1
            for raft in sorted(pool, key=lambda r:r.length, reverse=True):
                self.scaffold_map.append("{}\t{}\t{}\t{}\t{}\t{}\n".format(self.name, i+1, pool.pooltype, k, raft.name, raft.length))
                chrom_sbs += raft.manifest
                k += 1

                chrend = chrstart + raft.length - 1
                raft_orientation = '+'
                if raft.start > raft.end:
                    raft_orientation = '-'
                comment = 'anchored'
                if pool.pooltype == 'orient':
                    comment = 'unoriented'
                elif pool.pooltype == 'order':
                    comment = 'unordered'
                markers = ' '.join(['{:.3f}'.format(x) for x in sorted(set(raft.marker_chain))]) + ' cM'
                self.agp.append("{}\t{}\t{}\t{}\tD\t{}\t1\t{}\t{}\t# {}\t{}\n".format("chr" + self.name, chrstart, chrend, chrpart, raft.revised_name, raft.length, raft_orientation, comment, markers))
                chrpart += 1
                self.agp.append("{}\t{}\t{}\t{}\tN\t100\tfragment\tno\n".format("chr" + self.name, chrend+1, chrend+100, chrpart))
                chrom_sbs += [r.SummaryBlock(raft.revised_name + "_Gap", chrend+1, gd.Block(raft.revised_name + "_Gap", chrend+1, chrend+100, self.name, -1), 1)]
                chrpart += 1
                chrstart = chrend + 101
        
        # Delete final gap
        del self.agp[-1]
        del chrom_sbs[-1]
        
        chrommap = defaultdict(int)
        gaps = defaultdict(int)
        firstgap = 0
        curgap = 0
        curcm = -1
        for sb in chrom_sbs:
            if sb.cm == -1:
                if curcm == -1:
                    firstgap = sb.length
                else:
                    curgap += sb.length
            else:
                if curcm == sb.cm:
                    chrommap[curcm] += curgap
                elif curcm != -1:
                    gaps[curcm] = curgap

                chrommap[sb.cm] += sb.length
                curgap = 0
                curcm = sb.cm

        gaps[curcm] = curgap

        self.chromosome_map.append("{}\t{}\t{}\t{}\t{}\n".format(self.name, -1, 1, firstgap, firstgap))
        cmstart = firstgap + 1
        for cm in sorted(chrommap):
            cmend = cmstart + chrommap[cm] - 1
            self.chromosome_map.append("{}\t{:.3f}\t{}\t{}\t{}\n".format(self.name, cm, cmstart, cmend, chrommap[cm]))
            if gaps[cm] > 0:
                gapend = cmend+gaps[cm]
                self.chromosome_map.append("{}\t{}\t{}\t{}\t{}\n".format(self.name, -1, cmend+1, gapend, gaps[cm]))
                cmstart = gapend + 1
            else:
                cmstart = cmend + 1

    def remove_empty_pools(self):
        for i in reversed(range(len(self.pools))):
            self.pools[i].cleanup()
            if not self.pools[i].rafts:
                del self.pools[i]

    def connect(self, mergeclass):
        p = 0
        if type(mergeclass) is not list:
            mergeclass = [mergeclass]

        while p < len(self.pools)-1:
            if not self.pools[p]:
                del self.pools[p]
                continue

            q = p + 1
            while q < len(self.pools):
                for mc in mergeclass:
                    self.pools[p].assemble(self.pools[q], mc)
                q += 1
            p = self.split(p)
            p += 1

    def split(self, p):
        ordered_rafts = [raft for raft in self.pools[p] if raft.ordered]
        if not ordered_rafts or len(ordered_rafts) == len(self.pools[p]):
            return p

        self.pools.insert(p+1, pl.Pool(self))
        for raft in ordered_rafts:
            self.pools[p+1].add(raft)
            self.pools[p].remove(raft)
        return p+1

class Marker:
    def __init__(self, cm, prev_cm=-1, next_cm=-1):
        self.cm = cm
        self.prev_cm = prev_cm
        self.next_cm = next_cm

    def __repr__(self):
        return('{}-({},{})'.format(self.cm, self.prev_cm, self.next_cm))

    def update_previous(self, prev_cm):
        if prev_cm != -1:
            self.prev_cm = prev_cm
    
    def update_next(self, next_cm):
        if next_cm != -1:
            self.next_cm = next_cm


if __name__ == '__main__':
    
    pass