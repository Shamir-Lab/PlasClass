# Utility functions for the plasclass module

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
nt_bits = {'A':0,'C':1,'G':2,'T':3}

import numpy as np
import itertools

def readfq(fp): # this is a generator function
    ''' Adapted from https://github.com/lh3/readfq
    '''
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def get_rc(seq):
    ''' Return the reverse complement of seq
    '''
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])


def mer2bits(kmer):
    ''' convert kmer to bit representation
    '''
    bit_mer=nt_bits[kmer[0]]
    for c in kmer[1:]:
        bit_mer = (bit_mer << 2) | nt_bits[c]
    return bit_mer


def count_kmers(args_array):
    ''' Count the k-mers in the sequence
        Return a dictionary of counts
        Assumes ks is sorted
    '''
    ret_ind, seq, ks, kmer_inds, vec_lens, shared_list = args_array

    kmer_counts = {k:np.zeros(vec_lens[k]) for k in ks}

    k_masks = [2**(2*k)-1 for k in ks]
    ind=0
    bit_mers = [0 for k in ks]

    # get the first set of kmers
    while True:
        found = True
        for i,k in enumerate(ks):
            try:
                bit_mers[i] = mer2bits(seq[ind:ind+k])
                kmer_counts[k][kmer_inds[k][bit_mers[i]]] += 1.
            except:
                ind += 1
                found = False
                break
        if found == True:
            break

    # count all other kmers
    while ind<len(seq)-ks[-1]: # iterate through sequence until last k-mer for largest k
        for i,k in enumerate(ks):
            try:
                c = nt_bits[seq[ind+k]]
                bit_mers[i] = ((bit_mers[i]<<2)|c)&k_masks[i]
                kmer_counts[k][kmer_inds[k][bit_mers[i]]] += 1.
            except: # out of alphabet
                ind += 2 # pass it and move on to the next
                # get the next set of legal kmers
                while ind<=len(seq)-ks[-1]:
                    found = True
                    for i2,k2 in enumerate(ks):
                        try:
                            bit_mers[i2] = mer2bits(seq[ind:ind+k2])
                            kmer_counts[k2][kmer_inds[k2][bit_mers[i2]]] += 1.
                        except:
                            ind += 1
                            found = False
                            break
                    if found == True:
                        ind -= 1 # in next step increment ind
                        break
        ind += 1 # move on to next letter in sequence

    # count the last few kmers
    end = len(ks)-1
    for i in range(len(seq)-ks[-1]+1,len(seq)-ks[0]+1):
        for k in ks[:end]:
            kmer = seq[i:i+k]
            try:
                kmer_counts[k][kmer_inds[k][mer2bits(kmer)]] += 1.
            except:
                pass
        end -= 1

    #normalise counts
    kmer_freqs = np.zeros(sum([vec_lens[k] for k in ks]))
    ind = 0
    for k in ks:
        counts_sum = np.sum(kmer_counts[k])
        if counts_sum != 0:
            kmer_counts[k] = kmer_counts[k]/float(counts_sum)
        kmer_freqs[ind:ind+vec_lens[k]] = kmer_counts[k]
        ind += vec_lens[k]

    shared_list[ret_ind] = kmer_freqs


def compute_kmer_inds(ks):
    ''' Get the indeces of each canonical kmer in the kmer count vectors
    '''

    kmer_inds = {k: {} for k in ks}
    kmer_count_lens = {k: 0 for k in ks}

    alphabet = 'ACGT'
    for k in ks:
        all_kmers = [''.join(kmer) for kmer in itertools.product(alphabet,repeat=k)]
        all_kmers.sort()
        ind = 0
        for kmer in all_kmers:
            bit_mer = mer2bits(kmer)
            rc_bit_mer = mer2bits(get_rc(kmer))
            if rc_bit_mer in kmer_inds[k]:
                kmer_inds[k][bit_mer] = kmer_inds[k][rc_bit_mer]
            else:
                kmer_inds[k][bit_mer] = ind
                kmer_count_lens[k] += 1
                ind += 1
    return kmer_inds, kmer_count_lens
