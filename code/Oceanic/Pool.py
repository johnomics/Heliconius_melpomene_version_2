#!/usr/bin/env python3

from itertools import chain
from collections import defaultdict

class Pool:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.rafts = set()
        self.anchored = False # Used by NodeMerge to fix orient pool still at one marker only

    def __repr__(self):
        output = 'Chromosome {}\nType: {}\n'.format(self.chromosome.name, self.pooltype)
        for raft in self:
            if len(raft) == 0:
                continue
            output += repr(raft)
            output += '-----\n'
        output += '====='
        
        return output
    
    def __iter__(self):
        return iter(self.rafts)
    
    def __len__(self):
        return len(self.rafts)

    def add(self, raft):
        self.rafts.add(raft)

    def remove(self, raft):
        self.rafts.remove(raft)
    
    @property
    def scaffolds(self):
        return [scaffold for r in self.rafts for scaffold in r.scaffolds]
    
    @property
    def marker_chain(self):
        for raft in sorted(self.rafts, key=lambda x: x.length, reverse=True):
            return raft.marker_chain
        return ''
    
    @property
    def pooltype(self):
        nonempty_rafts = [r for r in self.rafts if len(r) > 0]
        if self.anchored:
            return 'ok'
        elif len(self.marker_chain) == 1:
            if len(nonempty_rafts) > 1:
                return 'order'
            else:
                return 'orient'
        else:
            if len(nonempty_rafts) > 1:
                return 'overlap'
            else:
                return 'ok'

    def cleanup(self):
        for raft in [raft for raft in self.rafts if not raft.logs]:
            self.rafts.remove(raft)

    def assemble(self, other, mergeclass, options=None):
        merger = mergeclass(self, other, options)
        for a in self.rafts:
            if not a:
                continue

            repeat = True
            seen = defaultdict(lambda:defaultdict(int))
            while repeat:
                repeat = False
                for b in other.rafts:
                    if not b or repr(a) == repr(b) or a in seen and b in seen[a]:
                        continue

                    seen[a][b] = 1
                    if merger.bridge(a, b):
                        del seen[a]
                        repeat = True
                        break
                if merger.merge():
                    del seen[a]
                    repeat = True
        other.cleanup()


    def extend(self):
        for raft in self.rafts:
            raft.extend()

    def write(self):
        scaffolds = []
        for raft in self:
            scaffolds.append(raft.write())
        return scaffolds

if __name__ == '__main__':
    pass
