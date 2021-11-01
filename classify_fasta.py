#!/usr/bin/env python
###
# Provide a command line script to classify sequences in a fasta file
###

from plasclass import plasclass_utils as utils
from plasclass import plasclass

import argparse
import sys

def parse_user_input():

    parser = argparse.ArgumentParser(
        description=
        'classify_fasta classifies the sequences in a fasta file as plasmid origin or not',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-f','--fasta',
     help='fasta file of the sequences to be classified',
     required=True, type=str
        )
    parser.add_argument('-o','--outfile',
     help='output file prefix',
     required=False, type=str
        )
    parser.add_argument('-p','--num_processes',
     help='Number of processes to use',
     required=False, type=int, default=8
        )

    return parser.parse_args()


def main(args):
    ''' Main function that classifies the sequences
    '''
    infile = args.fasta
    if args.outfile: outfile = args.outfile
    else: outfile = infile + '.probs.out'
    n_procs = args.num_processes

    c = plasclass.plasclass(n_procs)
    seq_names = []
    seqs = []
    print("Reading {} in batches of 100k sequences".format(infile))
    i = 0
    n_uppers = 0
    fp = open(infile)
    with open(outfile,'w') as o:
        for name, seq, _ in utils.readfq(fp):
            seq_names.append(name)
            if not seq.isupper():
                print("WARNING: sequence of {} converted to uppercase".format(name), file=sys.stderr)
                n_uppers += 1
                seq = seq.upper()
            seqs.append(seq)
            i += 1
            if i % 100000 == 0:
                print("Read {} sequences".format(i))
                probs = c.classify(seqs)
                for j,p in enumerate(probs):
                    o.write(seq_names[j] + '\t' + str(p) + '\n')
                seq_names = []
                seqs = []


        # last bunch of sequences:
        print("Read {} sequences".format(i))
        probs = c.classify(seqs)
        for j,p in enumerate(probs):
            o.write(seq_names[j] + '\t' + str(p) + '\n')
    fp.close()

    if n_uppers > 0:
        print('WARNING: {} sequences converted to uppercase, which is required for analysis'.format(n_uppers))
    print("Finished classifying")
    print("Class scores written in: {}".format(outfile))

if __name__=='__main__':
    args = parse_user_input()
    main(args)
