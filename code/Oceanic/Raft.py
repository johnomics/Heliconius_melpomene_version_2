#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Raft:
    def __init__(self, bl_id, scaffold, start, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.id = bl_id
        self.logs = []
        self.manifest = []
        self.append(scaffold, start, 1)
        self.bridges = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))

    def __repr__(self):
        output = ''
        output += self.scaffold + "\n"
        if not self.marker_chain.ok:
            output += 'To fix:\n'
        if self.scaffold in self.genome.offcuts:
            output += 'Offcut to {}\n'.format(','.join(self.genome.offcuts[self.scaffold]))
        if self.name in self.genome.haplotypes:
            output += 'Haplotype {}\n'.format(self.genome.haplotypes[self.name])
        output += '\n'.join([repr(m) for m in self.manifest]) + '\n'
        output += 'Length: {}\n'.format(self.length)
        return output

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.logs)
    
    def __len__(self):
        return len(self.logs)
        
    def __contains__(self, var):
        for log in self.logs:
            if log[0] == var or log[1] == var:
                return True
        return False

    @property
    def offcuts(self):
        if self.scaffold in self.genome.offcuts:
            return self.genome.offcuts[self.scaffold]
        return None

    @property
    def scaffolds(self):
        scaffolds = []
        for m in self.manifest:
            if m.scaffold not in scaffolds:
                scaffolds.append(m.scaffold)
        return scaffolds

    @property
    def scaffold(self):
        return '_'.join(self.scaffolds)
    
    def scaffold_range(self, scaffold):
        scaffold_sbs = [m for m in self.manifest if m.scaffold == scaffold]
        return scaffold_sbs[0].start, scaffold_sbs[-1].end

    @property
    def name(self):
        summarylist = []
        cur_scaffold = ""
        for m in self.manifest:
            if m.scaffold != cur_scaffold:
                cur_scaffold = m.scaffold
                summarylist.append([])
            summarylist[-1].append((m.scaffold, m.start, m.end))
        names = []
        for scaffold in summarylist:
            names.append('{}_{}_{}'.format(scaffold[0][0], scaffold[0][1], scaffold[-1][2]))
        name = '-'.join(names)
        return name

    @property
    def start(self):
        return self.manifest[0].start
    
    @property
    def last(self):
        return self.logs[-1][1]

    @property
    def end(self):
        return self.manifest[-1].end

    @property
    def length(self):
        length = 0
        for m in self.manifest:
            length += m.length
        return length

    @property
    def revised_name(self):
        if self.name not in self.genome.revised_names:
            self.genome.revised_names[self.name] = self.genome.revised + "{:02d}".format(int(self.chromosome.name)) + "{:03d}".format(self.chromosome.scaffold_count)
            self.chromosome.scaffold_count += 1
        return self.genome.revised_names[self.name]

    @property
    def sequence(self):
        sequence = sum([h.sequence for h in self.hooks], Seq("", generic_dna))
        sequence.id = self.revised_name
        sequence.description = "length={}".format(len(sequence))

        return sequence

    @property
    def marker_chain(self):
        mc = MarkerChain(self.manifest)
        return mc

    @property
    def ordered(self):
        return len(self.marker_chain) > 1

    @property
    def hooks(self):
        self._hooks = self.forge_hooks()
        return self._hooks

    def empty(self):
        self.logs = []
        self.manifest = []
        self.update()

    def discard(self, reason):
        self.genome.refuse.append((self.summary(), reason))
        self.empty()

    def replace(self, logs):
        self.empty()
        for scaffold, start, direction in logs:
            self.append(scaffold, start, direction)
        self.update()

    def reverse(self):
        newlogs = deepcopy(self.logs)
        self.empty()
        for scaffold, start, direction in reversed(newlogs):
            self.append(scaffold, start, direction*-1)

    def append(self, scaffold, start, direction=1):
        self.logs = self.logs + [(scaffold, start, direction)]
        self.manifest = self.manifest + [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)]
        self.update()

    def prepend(self, scaffold, start, direction=1):
        self.logs = [(scaffold, start, direction)] + self.logs
        self.manifest = [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)] + self.manifest
        self.update()

    def get_log_start(self, pos):
        for log in self.logs:
            if pos >= log[1] and pos <= self.genome.blocks[log[0]][log[1]].end:
                return log[1]
        else:
            return None

    def merge(self, other, before=False):
        if before:
            for scaffold, start, direction in reversed(other.logs):
                self.prepend(scaffold, start, direction)
        else:
            for scaffold, start, direction in other.logs:
                self.append(scaffold, start, direction)

    def update(self):
        self.remove_duplicates()
        self.collapse()
        self._hooks = self.forge_hooks()

    def remove_duplicates(self):
        dupes = defaultdict(int)
        new_logs = []
        for log in self.logs:
            dupes[log] += 1
            if dupes[log] <= 1:
                new_logs.append(log)
        self.logs = new_logs
        

        dupes = defaultdict(int)
        new_manifest = []
        for sb in self.manifest:
            sbkey = '{}_{}'.format(sb.scaffold, sb.start)
            dupes[sbkey] += 1
            if dupes[sbkey] <= 1:
                new_manifest.append(sb)
        self.manifest = new_manifest

        
    def collapse_consecutive(self):
        newsummary = []
        for i in range(0, len(self.manifest)):
            sbi = self.manifest[i]
            if len(newsummary) == 0:
                newsummary.append(sbi)
            else:
                last = newsummary[-1]
                if (last.scaffold==sbi.scaffold and
                    abs(last.end-sbi.start) == 1 and
                    last.cm == sbi.cm):
                    last.end = sbi.end
                else:
                    newsummary.append(sbi)

        self.manifest = newsummary

    def collapse_trios(self):
        if len(self.manifest) < 3:
            return

        newsummary = []
        i = 0
        while i < len(self.manifest):
            sbi = self.manifest[i]
            if i < len(self.manifest)-2:
                sbj = self.manifest[i+1]
                sbk = self.manifest[i+2]
                if (sbi.scaffold == sbj.scaffold and sbi.scaffold == sbk.scaffold and
                    abs(sbi.end - sbj.start) == 1 and abs(sbj.end - sbk.start) == 1 and
                    sbi.cm == sbk.cm             and sbj.cm == -1):
                        sbi.end = sbk.end
                        newsummary.append(sbi)
                        i += 3
                        continue

            newsummary.append(sbi)
            i += 1
        
        self.manifest = newsummary

    def collapse(self):
        self.collapse_consecutive()
        self.collapse_trios()

    def extend(self):
        self.extend_dir(0)
        self.extend_dir(-1)
        self.update()

    def extend_dir(self, item):
        first_scaffold, first_start, direction = self.logs[item]
        ext_start = first_start
        extend = 0
        starts_to_extend = []
        scaffold = self.genome.blocks[first_scaffold]
        chromosome = scaffold[first_start].chromosome

        while True:
            if item == 0 and direction == 1 or item == -1 and direction == -1:
                ext_start = scaffold[ext_start].prev_block
            else:
                ext_start = scaffold[ext_start].next_block
        
            # Extending to ends is OK
            if ext_start == 0:
                break
        
            # If we reach another cM, abandon extension
            if scaffold[ext_start].cm != -1 or scaffold[ext_start].chromosome != '0' and scaffold[ext_start].chromosome != chromosome:
                starts_to_extend = []
                break
        
            starts_to_extend.append(ext_start)
        
        if starts_to_extend:
            if item == 0:
                for start in starts_to_extend:
                    self.prepend(first_scaffold, start, direction)
            else:
                for start in starts_to_extend:
                    self.append(first_scaffold, start, direction)

    def write(self):
        if self.genome.revised_fasta:
            SeqIO.write(self.sequence, self.genome.revised_fasta, "fasta")
            for sb in self.manifest:
                self.chromosome.revised_db.execute("insert into scaffold_map values (?,?,?,?,?,?)",
                      [self.chromosome.name, '{:.3f}'.format(sb.cm), sb.scaffold, min(sb.start, sb.end), max(sb.start, sb.end), sb.length])
        
        return self.summary()
    
    def summary(self):
        scaffolds = []
        for scaffold, start, direction in self.logs:
            scaffolds.append(self.genome.blocks[scaffold][start])
        return scaffolds


    def forge_hooks(self):
        scaffold = start = end = None
        hooks = []
        for sb in self.manifest:
            if sb.scaffold != scaffold:
                if scaffold is not None:
                    hooks.append(Hook(self, scaffold, start, end))
                start = sb.start
                scaffold = sb.scaffold
            end = sb.end
        hooks.append(Hook(self, scaffold, start, end))

        return hooks


class MarkerChain:
    def __init__(self, manifest=None, chain=None):
        self.chain = []
        self.ok = True
        if manifest:
            for m in manifest:
                if m.cm != -1:
                    self.chain.append(m.cm)
                    if len(self.chain)>1 and self.chain[-2] == self.chain[-1]:
                        del self.chain[-1]
        elif chain:
            self.chain = chain

        if not self.check():
            self.ok = False

    def __getitem__(self, i):
        return self.chain[i]

    def __contains__(self, key):
        for cm in self.chain:
            if cm == key:
                return True
        else:
            return False

    def __repr__(self):
        return repr(self.chain)

    def __len__(self):
        if self.ok:
            return len(self.chain)
        else:
            return 0

    def __add__(self, other):
        new = self.chain + other.chain
        collapsed = [new[0]]
        for i in range(1,len(new)):
            if new[i] != collapsed[-1]:
                collapsed.append(new[i])
        return MarkerChain(chain=collapsed)

    def check(self):
        
        if len(self.chain) <= 1:
            return True
        
        if self.chain[0] < self.chain[1]:
            direction = 1
        else:
            direction = -1
        
        for i in range(1, len(self.chain)-1):
            if self.chain[i] < self.chain[i+1]:
                this_dir = 1
            else:
                this_dir = -1
            if direction != this_dir:
                return False
        
        return True



class SummaryBlock:
    def __init__(self, scaffold, start, block, direction):
        self.scaffold = scaffold
        if direction == 1:
            self.start = start
            self.end = block.end
        elif direction == -1:
            self.start = block.end
            self.end = start
        self.cm = block.cm

    @property
    def length(self):
        return max(self.start,self.end)-min(self.start,self.end)+1

        
    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))


class Hook():
    def __init__(self, raft, scaffold, start, end, knots=None):
        self.raft = raft
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.faces = []
        self.knots = {} if knots is None else knots

    def __repr__(self):
        return "{0:16s}:{1:7d}-{2:7d} {3}".format(self.scaffold, self.start, self.end, self.faces)

    def __lt__(self, other):
        if self.scaffold < other.scaffold:
            return self
        elif self.start < other.start:
            return self
        else:
            return other

    def __iter__(self):
        return iter(self.knots)

    @property
    def length(self):
        if self.start < self.end:
            return self.end - self.start + 1
        else:
            return self.start - self.end + 1
    
    @property
    def sequence(self):
        if 'Gap' in self.scaffold:
            seq = 'N' * self.length
        elif self.start < self.end:
            seq = self.raft.genome.sequences[self.scaffold][self.start-1:self.end]
        else:
            seq = self.raft.genome.sequences[self.scaffold][self.end-1:self.start].reverse_complement()
        return seq

    def order(self, start, end):
        if start < end:
            return start, end
        else:
            return end, start

if __name__ == '__main__':
    pass
