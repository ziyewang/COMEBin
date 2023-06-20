# from https://github.com/dparks1134/UniteM/blob/master/unitem/common.py

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import ast
from collections import defaultdict

import scripts.seq_io as seq_io
import biolib.seq_tk as seq_tk

from checkm.defaultValues import DefaultValues

from scripts.unitem_defaults import *


def calculateN50L50M50(seqs):
    """Calculate N50 and M50 statistics."""

    seq_lens = [len(seq) for seq in seqs]

    thresholdN50 = sum(seq_lens) / 2

    seq_lens.sort(reverse=True)

    test_sum = 0
    L50 = 0
    N50 = 0
    for seq_len in seq_lens:
        L50 += 1
        test_sum += seq_len
        if test_sum >= thresholdN50:
            N50 = seq_len
            break

    M50 = len(seqs) - L50
    if test_sum != thresholdN50:
        M50 += 1

    return N50, L50, M50


def parse_checkm_bin_stats(checkm_dir):
    """Read bin statistics from file."""
    bin_stats_file = os.path.join(checkm_dir, 'storage', DefaultValues.BIN_STATS_EXT_OUT)

    bin_stats = {}
    with open(bin_stats_file, 'r') as f:
        for line in f:
            line_split = line.split('\t')
            bin_stats[line_split[0]] = ast.literal_eval(line_split[1])

    return bin_stats


def parse_bin_stats(profile_dir):
    """Parse genomic and assembly statistics for bins."""

    binning_methods_dir = os.path.join(profile_dir, BINNING_METHOD_DIR)
    bin_stats = {}
    for bm in os.listdir(binning_methods_dir):
        checkm_dir = os.path.join(binning_methods_dir, bm, CHECKM_BAC_DIR)

        bin_stats[bm] = parse_checkm_bin_stats(checkm_dir)

    return bin_stats


def read_bins(bin_dirs):
    """Read sequences in bins."""

    bins = defaultdict(lambda: defaultdict(set))
    contigs = {}
    contigs_in_bins = defaultdict(lambda: {})
    for method_id, (bin_dir, bin_ext) in bin_dirs.items():
        for bf in os.listdir(bin_dir):
            if not bf.endswith(bin_ext):
                continue

            bin_id = bf[0:bf.rfind(bin_ext)]
            if bin_id[-1] == '.':
                bin_id = bin_id[0:-1]
            bf_path = os.path.join(bin_dir, bf)

            for seq_id, seq in seq_io.read_seq(bf_path):
                bins[method_id][bin_id].add(seq_id)
                contigs[seq_id] = seq
                contigs_in_bins[seq_id][method_id] = bin_id

            if len(bins[method_id][bin_id]) == 0:
                print('Bin %s from %s is empty.' % (bf, method_id))

    return bins, contigs, contigs_in_bins


def bin_gc(bin):
    """Calculate GC-content of bin."""

    a_count = 0
    c_count = 0
    g_count = 0
    t_count = 0
    for seq in bin.values():
        a, c, g, t = seq_tk.count_nt(seq)

        a_count += a
        c_count += c
        g_count += g
        t_count += t

    total = (a_count + c_count + g_count + t_count)
    return (g_count + c_count) / total
