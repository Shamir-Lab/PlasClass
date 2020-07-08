#!/usr/bin/env python

##
# Provide code to train new models
##

import argparse

import numpy as np
import os

import random

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from joblib import dump

import multiprocessing as mp
from multiprocessing import Manager

import plasclass_utils as utils

def parse_user_input():

    parser = argparse.ArgumentParser(
        description=
        'train.py trains models for the PlasClass classifier'
        )
    parser.add_argument('-p','--plasmid',
     help='plasmid file - file of plasmid reference sequences',
     required=True, type=str
     )
    parser.add_argument('-c','--chromosome',
     help='chromosome file - file of chromosome reference sequences',
     required=True, type=str
     )
    parser.add_argument('-o','--outdir',
     help='output directory',
     required=True, type=str
     )
    parser.add_argument('-n','--num_processes',
     help='number of processes to use',
     required=False, type=int, default=16
     )
    parser.add_argument('-k','--kmers',
     help='comma-separated list of k-mer lengths',
     required=False, type=str, default='3,4,5,6,7'
     )
    parser.add_argument('-l','--lengths',
     help='comma-separated list of sequence length bins',
     required=False, type=str, default='1000,10000,100000,500000'
     )

    return parser.parse_args()

def get_seq_lengths(infile):
    ''' Read in all the fasta entries,
        return arrays of the headers, and sequence lengths
    '''
    sequence_names = []
    sequence_lengths = []
    fp = open(infile)
    for name, seq, _ in utils.readfq(fp):
        sequence_names.append(name)
        sequence_lengths.append(len(seq))
    fp.close()

    return sequence_names, sequence_lengths


def get_num_frags(seq_lengths, length, coverage=5):
    ''' Compute how many fragments to generate
    '''
    # filter out sequences that are significantly shorter than the length
    filtered_seq_lengths = [l for l in seq_lengths if l > 0.85*length]
    tot_seq_len = sum(filtered_seq_lengths)

    num_frags_for_cov = int(np.ceil(tot_seq_len*coverage/float(length)))
    num_frags = min(90000,num_frags_for_cov)
    return num_frags


def get_start_inds(seq_names, seq_lengths, num_frags, length):
    ''' Randomly simulate fragments of a specific length from the sequences
    '''
    # filter out sequences that are significantly shorter than the length
    filtered_seq_names = [seq_names[i] for i,v in enumerate(seq_lengths) if v > 0.85*length]
    filtered_seq_lengths = [l for l in seq_lengths if l > 0.85*length]

    tot_seq_len = sum(filtered_seq_lengths)

    length_fractions = [float(l)/float(tot_seq_len) for l in filtered_seq_lengths]

    inds_dict = {}
    for name in filtered_seq_names:
        inds_dict[name] = []

    for i in range(num_frags):
        # choose genome
        seq_ind = np.random.choice(len(filtered_seq_names),p=length_fractions)
        seq_len = filtered_seq_lengths[seq_ind]
        seq_name = filtered_seq_names[seq_ind]
        # choose start index in the genome
        if seq_len < length: # just take the whole thing
            inds_dict[seq_name].append(0)
        else:
            start_ind = random.randint(0,seq_len - 1)
            inds_dict[seq_name].append(start_ind)

    return inds_dict


def get_seqs(infile, inds_dict, l):
    ''' Create array of the sequences
    '''
    seqs = []
    fp = open(infile)
    for name,seq,_ in utils.readfq(fp):
        start_inds = inds_dict.get(name, [])
        for start_ind in start_inds:
            frag = seq[start_ind:start_ind+l]
            if len(frag) < l and len(seq) > l:
                frag += seq[:l-len(frag)]
            seqs.append(frag)
    fp.close()
    return seqs


def train(plasfile, chromfile, outdir, num_procs, ks=[3,4,5,6,7], lens=[1000,10000,100000,500000]):
    ''' Train PlasClass models
    '''
    print("Starting PlasClass training")
    print("Getting reference lengths")
    chrom_names, chrom_lengths = get_seq_lengths(chromfile)
    plas_names, plas_lengths = get_seq_lengths(plasfile)
    for l in lens:
        coverage=5 # TODO: make this command line option
        num_frags = get_num_frags(plas_lengths,l,coverage)

        print("Sampling {} fragments for length {}".format(num_frags,l))
        plas_start_inds = get_start_inds(plas_names, plas_lengths, num_frags, l)
        chrom_start_inds = get_start_inds(chrom_names, chrom_lengths, num_frags, l)
        plas_seqs = get_seqs(plasfile, plas_start_inds, l)
        chrom_seqs = get_seqs(chromfile, chrom_start_inds, l)

        print("Getting k-mer frequencies")
        kmer_inds, kmer_count_lens = utils.compute_kmer_inds(ks)

        pool = mp.Pool(num_procs)
        plas_list=Manager().list()
        for cur in np.arange(len(plas_seqs)):
            plas_list.append(0)
        pool.map(utils.count_kmers, [[ind,s, ks, kmer_inds, kmer_count_lens, plas_list] \
                                        for ind,s in enumerate(plas_seqs)])
        plas_freqs = np.array(plas_list)

        chrom_list=Manager().list()
        for cur in np.arange(len(chrom_seqs)):
            chrom_list.append(0)
        pool.map(utils.count_kmers, [[ind, s, ks, kmer_inds, kmer_count_lens, chrom_list] \
                                        for ind,s in enumerate(chrom_seqs)])
        chrom_freqs = np.array(chrom_list)

        pool.close()

        print("Learning classifier")
        plas_labels = np.ones(plas_freqs.shape[0])
        chrom_labels = np.zeros(chrom_freqs.shape[0])
        data = np.concatenate((plas_freqs,chrom_freqs))
        labels = np.concatenate((plas_labels, chrom_labels))
        scaler = StandardScaler().fit(data)
        scaled = scaler.transform(data)
        clf = LogisticRegression(solver='liblinear').fit(scaled,labels)

        print("Saving classifier")
        clf_name = 'm'+str(l)
        scaler_name = 's'+str(l)
        dump(clf, os.path.join(outdir,clf_name))
        dump(scaler, os.path.join(outdir,scaler_name))


def main(args):
    ''' Run the training code
    '''
    plasfile = args.plasmid
    chromfile = args.chromosome
    ks = [int(k) for k in args.kmers.split(',')]
    lens = [int(s) for s in args.lengths.split(',')]
    num_procs = args.num_processes
    outdir = args.outdir

    train(plasfile,chromfile,outdir,num_procs,ks,lens)

if __name__=='__main__':
    args = parse_user_input()
    main(args)
