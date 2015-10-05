#!/usr/bin/python3 -u

import os
import sys
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description='''Filter within-genome hits from chain file (gzipped)
    
        -c chainfile (gzipped)
        -p prefix
        -o output
        -l length_threshold
        ''')

    parser.add_argument('-c', '--chainfile', type=str, required=True)
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-l', '--length', type=int, required=False)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    outgz = gzip.open(args.output, 'wb')
    with gzip.open(args.chainfile, 'rt') as c:
        valid = True
        for line in c:
            if line.startswith('#'):
                continue
            if line.startswith('chain'):
                chainargs = line.rstrip().split(' ')
                valid = not chainargs[2].startswith(args.prefix) and chainargs[7].startswith(args.prefix) or chainargs[2].startswith(args.prefix) and not chainargs[7].startswith(args.prefix)
                if valid and args.length:
                    valid = args.length < (int(chainargs[6])-int(chainargs[5])+1) and args.length < (int(chainargs[11])-int(chainargs[10])+1)
            if valid:
                outgz.write(bytes(line,'UTF-8'))