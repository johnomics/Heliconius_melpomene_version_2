#!/usr/bin/env python3

from copy import deepcopy
from collections import defaultdict, namedtuple
from math import sqrt
from operator import itemgetter

from . import Raft as r
from . import GenomeData as gd

class Merger():
    
    def order(self, start, end):
        if start < end:
            return start, end
        else:
            return end, start
    
    def get_pool_positions(self, this, other):
        i = j = None
        for n, p in enumerate(this.chromosome.pools):
            if this is p:
                i = n
            if other is p:
                j = n
        return i,j
    
    def do_join(self, a, b, bridge):
        a.replace(bridge)
        b.empty()


class MarkerMerge(Merger):
    
    def __init__(self, this, other, options=None):
        self.joins = []

    def merge(self):
        return False

    def bridge(self, a, b):
        if a.scaffold != b.scaffold:
            return None
        bridge = []
        scaffold = a.genome.blocks[a.scaffold]
    
        start, target, a_forward, b_forward = self.orient(a,b)
        this_cm = scaffold[start].cm
        direction = 1 if start < target else -1
        
        cm_dir = 0

        while True:
            if direction == 1:
                start = scaffold[start].next_block
            else:
                start = scaffold[start].prev_block

            if start == 0:
                break
            
            block_i_chr = scaffold[start].chromosome
            block_i_cm = scaffold[start].cm
            if block_i_chr == '0' or block_i_cm == -1 and block_i_chr == a.chromosome.name:
                bridge.append((a.scaffold, start, direction))
                continue
        
            if block_i_chr != a.chromosome.name: # This check comes after the cM check because chr==0 is OK
                break
        
            merge = False
            if block_i_cm == this_cm:
                merge = True
        
            if block_i_cm == a.chromosome.markers[this_cm].next_cm and cm_dir in (0, 1):
                this_cm = a.chromosome.markers[this_cm].next_cm
                cm_dir = 1
                merge = True
        
            
            if block_i_cm == a.chromosome.markers[this_cm].prev_cm and cm_dir in (0, -1):
                this_cm = a.chromosome.markers[this_cm].prev_cm
                cm_dir = -1
                merge = True
                
            
            if merge:
                if start == target:
                    bridge = self.merge_logs(bridge, a.logs)
                    bridge = self.merge_logs(bridge, b.logs)
                    self.do_join(a, b, bridge)
                    return True
                else:
                    bridge.append((a.scaffold, start, direction))
                
            if bridge and bridge[-1][1] == start:
                continue
    
            break
    
        return None

    def orient(self, a, b):
        edge_blocks = [a.start, a.last, b.start, b.last]
        if a.start < b.start:
            a_far, a_near, b_near, b_far = sorted(edge_blocks)
        else:
            b_far, b_near, a_near, a_far = sorted(edge_blocks)
        
        a_forward = a.start <= a.last
        b_forward = b.start <= b.last
        return a.get_log_start(a_near), b.get_log_start(b_near), a_forward, b_forward
    
    def merge_logs(self, this, other):
        if not this:
            return deepcopy(other)
    
        new = deepcopy(this) + deepcopy(other)
        new.sort(key=lambda x: x[1]) # Sort by start position
    
        this_dir = this[0][2]
        new = [(scaffold, start, this_dir) for scaffold, start, direction in new]
        if this_dir == -1:
            new.reverse()
    
        return new


class NodeMerge(Merger):
    def __init__(self, this, other, options=None):
        self.this = this
        self.other = other
        self.genome = this.genome
        self.connections = defaultdict(lambda: defaultdict(NodeMerge.Connection))
    
    Connection = namedtuple('Connection', 'node, anchor_raft, floating_raft, floating_pool, float_after')
    
    def bridge(self, a, b):
        i, j = self.get_pool_positions(self.this, self.other)
        if self.this.pooltype == 'ok' and self.other.pooltype == 'ok': # Will be handled by OKMerge
            return False

        if j not in [i, i+1]:   # Only consider same pool or neighbouring pools
            return False

        if self.this.pooltype == 'ok':
            anchor_pool, floating_pool = self.this, self.other
            anchor_raft, floating_raft = a, b
            float_after = True
        elif self.other.pooltype == 'ok':
            anchor_pool, floating_pool = self.other, self.this
            anchor_raft, floating_raft = b, a
            float_after = False
        else:
            return False    # Only process pairs with one OK pool

        anchor_scaffold = anchor_raft.scaffolds[0]
        if float_after:
            anchor_scaffold = anchor_raft.scaffolds[-1]

        for floating_scaffold in floating_raft.scaffolds:
            if floating_scaffold in self.connections and anchor_scaffold in self.connections[floating_scaffold]:
                continue
            directions = {}
            for node in sorted(self.genome.nodes[anchor_scaffold][floating_scaffold], key = lambda x:x.tlen, reverse=True): # Longest hit first
                scaffold_start, scaffold_end = anchor_raft.scaffold_range(anchor_scaffold)
                node.set_status(scaffold_start, scaffold_end, floating_raft.start, floating_raft.end)
                if 'connect' in node.status:
                    self.connections[anchor_scaffold][floating_scaffold] = NodeMerge.Connection(node, anchor_raft, floating_raft, floating_pool, float_after)
                    break

        return False

    def merge(self):
        if len(self.connections) == 0:
            return False

        for anchor_scaffold in self.connections:

            if len(self.connections[anchor_scaffold]) > 1:  # Only merge when only one possible connecting scaffold
                continue

            floating_scaffold = list(self.connections[anchor_scaffold])[0]
            conn = self.connections[anchor_scaffold][floating_scaffold]
            reverse_float = conn.anchor_raft.start < conn.anchor_raft.end and conn.node.direction == -1 or conn.anchor_raft.start > conn.anchor_raft.end and conn.node.direction == 1
            if conn.floating_pool.pooltype == 'orient': # Reverse if necessary and anchor; OKMerge will take care of the merging
                if reverse_float:
                    conn.floating_raft.reverse()
                conn.floating_pool.anchored = True
                return False
            
            elif conn.floating_pool.pooltype == 'order':
                anchor_tip = conn.anchor_raft.end if conn.float_after else conn.anchor_raft.start
                if 'connect_a' in conn.node.status and ('af' in conn.node.status and anchor_tip < conn.node.tend or 'ar' in conn.node.status and anchor_tip > conn.node.tstart):
                    return False
                if 'connect_b' in conn.node.status and ('af' in conn.node.status and anchor_tip > conn.node.tstart or 'ar' in conn.node.status and anchor_tip < conn.node.tend):
                    return False

                if reverse_float:
                    conn.floating_raft.reverse()
                newgap = self.genome.add_gap()
                if conn.float_after:
                    conn.anchor_raft.append(newgap.scaffold, newgap.start, 1)
                    conn.anchor_raft.merge(conn.floating_raft)
                else:
                    conn.anchor_raft.prepend(newgap.scaffold, newgap.start, 1)
                    conn.anchor_raft.merge(conn.floating_raft, before=True)
                conn.floating_raft.empty()

                del self.connections[anchor_scaffold][floating_scaffold]
                if len(self.connections[anchor_scaffold]) == 0:
                    del self.connections[anchor_scaffold]
                return True
        
        return False


class OKMerge(Merger):
    def __init__(self, this, other, options=None):
        self.this = this
        self.other = other
        self.genome = this.genome
        
    def bridge(self, a, b):
        if self.this.pooltype != 'ok' or self.other.pooltype != 'ok':
            return False

        i, j = self.get_pool_positions(self.this, self.other)

        for n in range(i+1, j):
            if len(self.this.chromosome.pools[n]) != 0 or self.this.chromosome.pools[n].pooltype in ['other', 'orient']:
                return False
        
        newgap = self.genome.add_gap()
        a.append(newgap.scaffold, newgap.start, 1)
        
        a.merge(b)
        b.empty()

        return True
        
    def merge(self):
        return False


if __name__ == '__main__':
    pass