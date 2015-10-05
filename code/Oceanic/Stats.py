#!/usr/bin/env python3

class PoolStat:
    def __init__(self, count=0, length=0):
        self.count = count
        self.length = length
    
    def __add__(self, other):
        self.count += other.count
        self.length += other.length
        return self

class Stats:
    def __init__(self, name):
        self.name = name
        self.pool_num = 0
        self.log_num = 0
        self.length = 0
        self.pool_stats = {'orient':PoolStat(), 'order':PoolStat(), 'overlap':PoolStat(), 'error':PoolStat(), 'ok':PoolStat()}
        self.log_count = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}
        self.raft_count = {'orient':0, 'order':0, 'overlap':0, 'error':0, 'ok':0}

    def __add__(self, other):
        self.pool_num += other.pool_num
        self.log_num += other.log_num
        self.length += other.length
        for gt in self.pool_stats:
            self.pool_stats[gt] += other.pool_stats[gt]
            self.log_count[gt] += other.raft_count[gt]
            self.raft_count[gt] += other.raft_count[gt]
        return self

    def __repr__(self):
        output = "{}: {:6d} blocks in {:4d} groups, {:9d} bp\n".format(self.name, self.log_num, self.pool_num, self.length)
        output += "Type\tGroups\tLists\tBlocks\tLength\n"
        for gt in 'ok', 'orient', 'order', 'overlap', 'error':
            output += "{}\t{:6d}\t{:6d}\t{:6d}\t{:10d}\n".format(
                  gt, self.pool_stats[gt].count, self.raft_count[gt],
                  self.log_count[gt], self.pool_stats[gt].length)
        return output


def genome(lengths):
    scaffolds=len(lengths)
    genome_length = sum(lengths)
    lengths.sort(reverse=True)
    N50_threshold = genome_length / 2
    
    N50_sum = 0
    N50 = 0
    for length in lengths:
        N50_sum += length
        if N50_threshold < N50_sum:
            N50 = length
            break
    
    return "Scaffolds: {}\tLength: {}\tN50: {}".format(scaffolds, genome_length, N50)


if __name__ == '__main__':
    pass