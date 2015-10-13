#!/usr/bin/env python3

# John Davey jd626@cam.ac.uk

import argparse
import re
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''Generate AGP file from one-line FASTA file.
    
    -f FASTA file''')

parser.add_argument('-f', '--fastafile', type=str, required=True)

args=parser.parse_args()

def make_part(seq, scaffold, pos, partnum, type):
    length = len(seq)
    end = pos + length - 1
    part = '\t'.join([scaffold, str(pos), str(end), str(partnum), type])

    if type is 'W':
        contig_name = scaffold+'.'+str(partnum)
        part = '\t'.join([part, contig_name, str(1), str(length), '+\n'])
    if type is 'N':
        part = '\t'.join([part, str(length), 'fragment', 'yes\n'])

    pos += len(seq)
    partnum += 1
    return(part, pos, partnum)


header_re = re.compile(r'^>(.+)$')
seq_re = re.compile(r'([ACGT]+)?(N+)?', re.IGNORECASE)

try:
    agpname = args.fastafile + ".agp"
    with open(agpname, 'w') as agp:
        for record in SeqIO.parse(args.fastafile, "fasta"):
            scaffold = record.id
            pos = 1
            partnum = 1
            for m in seq_re.finditer(str(record.seq)):
                (seq, ns) = m.groups()
                if seq:
                    (part, pos, partnum) = make_part(seq, scaffold, pos, partnum, 'W')
                    agp.write(part)
                if ns:
                    (part, pos, partnum) = make_part(ns, scaffold, pos, partnum, 'N')
                    agp.write(part)

except IOError:
        print("Cannot open file " + args.fastafile + "!")

    