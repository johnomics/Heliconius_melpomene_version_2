#!/usr/bin/python3 -u

import os
import sys
import argparse
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

import Oceanic.GenomeData as gd
import Oceanic.Chromosome as chrom
import Oceanic.Stats as stats

from Bio import SeqIO
from Bio.Seq import Seq

def get_genome_stats(chromosomes):
    genome_stats = stats.Stats("Genome")

    for chromosome in chromosomes:
        genome_stats += chromosomes[chromosome].stats

    print(genome_stats)

def load_map(genome):

    print("Loading map...")

    mapped_blocks = 0
    mapped_blocks_length = 0
    placed_blocks = 0
    placed_blocks_length = 0

    genome.db.execute("select distinct chromosome from scaffold_map order by chromosome")
    chromosomes = {chromosome_name: chrom.Chromosome(chromosome_name, genome) for chromosome_name, in genome.db.fetchall() if chromosome_name != 0}

    for name in chromosomes:
        mapped_blocks += chromosomes[name].mapped_blocks
        mapped_blocks_length += chromosomes[name].mapped_blocks_length
        placed_blocks += chromosomes[name].placed_blocks
        placed_blocks_length += chromosomes[name].placed_blocks_length

    print("Map blocks {} length {} of which {} blocks placed, length {}".format(mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length))
    return chromosomes

def collapse_unmapped_blocks(unmapped, stats):
    outblocks = []
    prevend = -1
    curstart = 0
    
    for start in unmapped:
        stats['blocks'] += 1
        stats['length'] += unmapped[start].length

    i = 0
    while i < len(unmapped)-1:
        starts = sorted(unmapped)
        ui = unmapped[starts[i]]
        ui.prev_block, ui.next_block = 0, 0
        uj = unmapped[starts[i+1]]
        uj.prev_block, uj.next_block = 0, 0
        if ui.chromosome == uj.chromosome and ui.cm == uj.cm and ui.end+1 == uj.start:
            ui.end = uj.end
            del unmapped[starts[i+1]]
            continue
        i += 1


def write_unmapped_scaffold(outblocks, scaffold, stats, unmapped_output, genome, chromosomes):

    seq = Seq('')
    dbblocks = []
    mapped = False
    scaffold_chromosomes = []
    for block in outblocks:
        if genome.revised_db:
            dbblocks.append([block.chromosome, block.cm, scaffold, block.start, block.end, block.length])                
        if genome.revised_fasta:
            seq += genome.sequences[scaffold][block.start-1:block.end]
        if block.chromosome != '0':
            scaffold_chromosomes.append(int(block.chromosome))


    stats['scaffolds'] += 1

    scaffold_chromosomes = set(scaffold_chromosomes)
    scaffold_start, scaffold_end = outblocks[0].start, outblocks[-1].end
    scaffold_length = scaffold_end - scaffold_start + 1
    stats['scaffold_length'] += scaffold_length

    if seq.features or len(scaffold_chromosomes) > 0 and len(dbblocks) > 1:
        unmapped_output.append([gd.Block(scaffold, scaffold_start, scaffold_end)])

        scaffold_name = genome.revised + "{:05d}".format(genome.revised_count)
        genome.revised_count += 1
        genome.revised_names["{}_{}_{}".format(scaffold, scaffold_start, scaffold_end)] = scaffold_name

        stats['written_scaffolds'] += 1
        stats['written_length'] += scaffold_length
        if genome.revised_db:
            for block in dbblocks:
                genome.revised_db.execute("insert into scaffold_map values (?,?,?,?,?,?)", block)
        if genome.revised_fasta:
            seq.description = "length={}".format(len(seq))
            seq.id = scaffold_name
            SeqIO.write(seq, genome.revised_fasta, "fasta")

        if len(scaffold_chromosomes) > 0:
            chrom = next(iter(scaffold_chromosomes))
            chr_unmapped_end = chromosomes[chrom].unmapped_start + scaffold_length - 1 
            chromosomes[chrom].agp.append("{}\t{}\t{}\t{}\tD\t{}\t1\t{}\t+\n".format("chr{}_unmapped".format(chrom), chromosomes[chrom].unmapped_start, chr_unmapped_end, chromosomes[chrom].unmapped_part, scaffold_name, scaffold_length))
            chromosomes[chrom].unmapped_part += 1
            chromosomes[chrom].agp.append("{}\t{}\t{}\t{}\tN\t100\tfragment\tno\n".format("chr{}_unmapped".format(chrom), chr_unmapped_end+1, chr_unmapped_end+100, chromosomes[chrom].unmapped_part))
            chromosomes[chrom].unmapped_part += 1
            chromosomes[chrom].unmapped_start = chr_unmapped_end + 101
    else:
        stats['discard_scaffolds'] += 1
        stats['discard_length'] += scaffold_length
        for dbblock in dbblocks:
            dblength = dbblock[5]
            partslength = 0
            for newpart in genome.newparts[scaffold]:
                for origpart in genome.origparts[newpart.oldname]:
                    if dbblock[2] == origpart.newname and (dbblock[3] <= origpart.newstart <= dbblock[4] or dbblock[3] <= origpart.newend <= dbblock[4]):
                        if origpart.parttype in ['active', 'retained']:
                            partslength += origpart.newend - origpart.newstart + 1
                            origpart.parttype = 'removed'
            if dblength != partslength:
                print(scaffold, dblength, partslength, dbblock)


def process_unmapped_scaffold(unmapped, scaffold, stats, genome, chromosomes):

    collapse_unmapped_blocks(unmapped, stats)
    
    unmapped_output = []
    outblocks = []
    starts = sorted(unmapped)
    for i, start in enumerate(starts):
        end = unmapped[start].end

        outblocks.append(unmapped[start])

        if i == len(starts)-1 or end+1 != unmapped[starts[i+1]].start:
            write_unmapped_scaffold(outblocks, scaffold, stats, unmapped_output, genome, chromosomes)
            outblocks = []

    return unmapped_output

def reassemble(chromosomes, genome, args):

    print("Reassembly...")

    pool = ThreadPool(args.threads)
    pool.map(lambda x: chromosomes[x].assemble(args), chromosomes.keys())
    get_genome_stats(chromosomes)
    
    gap_blocks = gap_length = total_blocks = total_length = 0
    for scaffold in genome.blocks:
        for start in genome.blocks[scaffold]:
            total_blocks += 1
            total_length += genome.blocks[scaffold][start].length
            if "Gap" in scaffold:
                gap_blocks += 1
                gap_length += genome.blocks[scaffold][start].length
    print("Added {} gaps, length {}".format(gap_blocks, gap_length))
    print("Total blocks now {}, length {}".format(total_blocks, total_length))

    assembly = []
    for chromosome in chromosomes:
        assembly += chromosomes[chromosome].write()


    assembled_blocks = assembled_length = 0
    for blocklist in assembly:
        for block in blocklist:
            assembled_blocks += 1
            assembled_length += block.length
            del genome.blocks[block.scaffold][block.start]
            if not genome.blocks[block.scaffold]:
                del genome.blocks[block.scaffold]

    print("Assembled {} blocks, length {}".format(assembled_blocks, assembled_length))

    remainder_blocks = remainder_length = 0
    for scaffold in genome.blocks:
        for start in genome.blocks[scaffold]:
            remainder_blocks += 1
            remainder_length += genome.blocks[scaffold][start].length
    print("Remainder: {} blocks, length {}".format(remainder_blocks, remainder_length))

    deleted_blocks = deleted_length = 0
    for (blocklist, reason) in genome.refuse:
        for block in blocklist:
            deleted_blocks += 1
            deleted_length += block.length
            for part in genome.newparts[block.scaffold]:
                if part.parttype not in  ['haplotype', reason] and \
                   (block.start >= part.newstart and block.end <= part.newend):
                    part.parttype = reason

            del genome.blocks[block.scaffold][block.start]
            if not genome.blocks[block.scaffold]:
                del genome.blocks[block.scaffold]

    print("Deleted {} blocks, length {}".format(deleted_blocks, deleted_length))


    unmapped_stats = {}
    unmapped_stats['blocks'] = 0
    unmapped_stats['length'] = 0
    unmapped_stats['scaffolds'] = 0
    unmapped_stats['scaffold_length'] = 0
    unmapped_stats['written_scaffolds'] = 0
    unmapped_stats['written_length'] = 0
    unmapped_stats['discard_scaffolds'] = 0
    unmapped_stats['discard_length'] = 0

    left_blocks = left_length = offcut_blocks = offcut_length = 0
    for scaffold in sorted(genome.blocks):
        if scaffold in genome.offcuts:
            for start in genome.blocks[scaffold]:
                offcut_blocks += 1
                offcut_length += genome.blocks[scaffold][start].length
                for part in genome.newparts[scaffold]:
                    for oldpart in genome.origparts[part.oldname]:
                        if part == oldpart:
                            oldpart.parttype = "offcut_unmapped_removed"
        else:
            length = 0
            unmapped_output = process_unmapped_scaffold(genome.blocks[scaffold], scaffold, unmapped_stats, genome, chromosomes)
            for unmapped_scaffold in unmapped_output:
                assembly.append(unmapped_scaffold)

    if genome.revised_db:
        genome.revised_conn.commit()

    if genome.revised_orig_tsv:
        for origscaffold in sorted(genome.origparts):
            for part in genome.origparts[origscaffold]:
                genome.revised_orig_tsv.write('{}\n'.format(repr(part)))

    if genome.revised_names and genome.revised_tsv:
        for name in sorted(genome.revised_names):
            parts = name.split('-')
            newstart = 1
            for p in parts:
                pscaffold, pstart, pend = p.split('_')
                pstart, pend = int(pstart), int(pend)
                if "Gap" in p:
                    newstart += pend - pstart + 1
                    continue
                if pstart < pend:
                    strand = 1
                    outstart, outend = pstart, pend
                else:
                    strand = -1
                    outstart, outend = pend, pstart
                length = outend - outstart + 1
                newend = newstart + length - 1
                output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\tactive\n".format(pscaffold, outstart, outend, genome.revised_names[name], newstart, newend, strand)
                genome.revised_tsv.write(output)
                newstart = newend + 1


    print("Removed {} unmapped offcuts, length {}".format(offcut_blocks, offcut_length))
    if unmapped_stats['scaffold_length'] != unmapped_stats['length']:
        print("unmapped stats don't match! {} {}".format(unmapped_stats['scaffold_length'], unmapped_stats['length']))
    print("Left over    {} blocks in {} scaffolds, length {}".format(unmapped_stats['blocks'], unmapped_stats['scaffolds'], unmapped_stats['length']))
    print("    of which {} scaffolds were written,   length {}".format(unmapped_stats['written_scaffolds'], unmapped_stats['written_length']))
    print("    and      {} scaffolds were discarded, length {}".format(unmapped_stats['discard_scaffolds'], unmapped_stats['discard_length']))

    return assembly


def output(assembly, chromosomes, revised):
    
    print("Output...")
    
    genome_length = 0
    scaffolds = 0
    scaffold_lengths = []
    blocks = 0
    for blockgroup in assembly:
        scaffold_length = 0
        for block in blockgroup:
            blocks += 1
            scaffold_length += block.length
        scaffolds += 1
        genome_length += scaffold_length
        scaffold_lengths.append(scaffold_length)
    print("Blocks:{}\t".format(blocks), stats.genome(scaffold_lengths))
    
    sm = open(revised+'_scaffold_map.tsv', 'w')
    sm.write("Chromosome\tPool\tType\tID\tScaffold\tLength\n")

    cm = open(revised+'_chromosome_map.tsv', 'w')
    cm.write("Chromosome\tcM\tStart\tEnd\tLength\n")

    agp = open(revised+'_chromosomes.agp', 'w')
    for c in chromosomes:
        for mapline in chromosomes[c].scaffold_map:
            chromosome, pool, pooltype, poolid, scaffold, length = mapline.split('\t')
            scaffold = genome.revised_names[scaffold]
            sm.write('\t'.join([chromosome, pool, pooltype, poolid, scaffold, length]))
        for mapline in chromosomes[c].chromosome_map:
            cm.write(mapline)

        if 'unmapped' in chromosomes[c].agp[-1] and '\tN\t' in chromosomes[c].agp[-1]:
            del chromosomes[c].agp[-1] # Final gap in unmapped scaffolds (Chromosome removes it from mapped scaffolds)
        for agpline in chromosomes[c].agp:
            agp.write(agpline)
    sm.close()
    cm.close()
    agp.close()


def get_args():
    parser = argparse.ArgumentParser(description='''Output new FASTA file based on linkage map.
    
        -d database
        -f FASTA file
        -g GFF file
        -e errors file
        -t threads
        -r revised
        -a haplomerger
        -o original TSV
        -p prefix
        -n nodes''')

    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-f', '--fasta', type=str, required=True)
    parser.add_argument('-g', '--gff', type=str, required=True)
    parser.add_argument('-e', '--errors', type=str)
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('-r', '--revised', type=str)
    parser.add_argument('-a', '--haplomerger', type=str)
    parser.add_argument('-o', '--original', type=str)
    parser.add_argument('-p', '--prefix', type=str)
    parser.add_argument('-n', '--nodes', type=str)

    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    genome = gd.GenomeData(args)

    chromosomes = load_map(genome)

    assembly = reassemble(chromosomes, genome, args)

    output(assembly, chromosomes, args.revised)