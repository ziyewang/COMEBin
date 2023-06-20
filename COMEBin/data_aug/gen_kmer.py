#!/usr/bin/env python
from __future__ import print_function
import os
import numpy as np
import pandas as pd
from itertools import product
from Bio import SeqIO
# optimized sliding window function from
# http://stackoverflow.com/a/7636587
from itertools import tee
from collections import Counter, OrderedDict
import pandas as p

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in range(i):
            next(el, None)
    return zip(*els)

def generate_feature_mapping(kmer_len):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash = {}
    counter = 0
    for kmer in product("ATGC",repeat=kmer_len):
        if kmer not in kmer_hash:
            kmer_hash[kmer] = counter
            rev_compl = tuple([BASE_COMPLEMENT[x] for x in reversed(kmer)])
            kmer_hash[rev_compl] = counter
            counter += 1
    return kmer_hash,counter

def generate_features_from_fasta(fasta_file,length_threshold,kmer_len,outfile):
    kmer_dict,nr_features = generate_feature_mapping(kmer_len)

    # Store composition vectors in a dictionary before creating dataframe
    composition_d = OrderedDict()
    contig_lengths = OrderedDict()
    for seq in SeqIO.parse(fasta_file,"fasta"):
        seq_len = len(seq)
        if seq_len <= length_threshold:
            continue
        contig_lengths[seq.id] = seq_len
        # Create a list containing all kmers, translated to integers
        kmers = [
            kmer_dict[kmer_tuple]
            for kmer_tuple
            in window(str(seq.seq).upper(), kmer_len)
            if kmer_tuple in kmer_dict
        ]
        kmers.append(nr_features-1)
        composition_v = np.bincount(np.array(kmers,dtype=np.int64))
        composition_v[-1]-=1
        composition_d[seq.id] = composition_v 
    df = p.DataFrame.from_dict(composition_d, orient='index', dtype=float)
    df.to_csv(outfile)

def run_gen_kmer(fasta_file, length_threshold, kmer_len):
    outfile = os.path.join(os.path.dirname(fasta_file), 'kmer_' + str(kmer_len) + '_f' + str(length_threshold) + '.csv')
    generate_features_from_fasta(fasta_file,length_threshold,kmer_len,outfile)

