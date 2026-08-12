from itertools import product

import numpy as np


_BASE_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
_BASE_CODE = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
_ASCII_TO_BASE = np.full(256, 255, dtype=np.uint8)
_ASCII_TO_BASE[[ord('A'), ord('a')]] = 0
_ASCII_TO_BASE[[ord('T'), ord('t')]] = 1
_ASCII_TO_BASE[[ord('G'), ord('g')]] = 2
_ASCII_TO_BASE[[ord('C'), ord('c')]] = 3
_ASCII_TO_BASE.setflags(write=False)


def _encode_kmer(kmer):
    encoded = 0
    for base in kmer:
        encoded = encoded * 4 + _BASE_CODE[base]
    return encoded


def build_feature_lookup(kmer_len=4):
    if kmer_len != 4:
        raise ValueError('COMEBin currently requires kmer_len=4.')

    lookup = np.full(4 ** kmer_len, 255, dtype=np.uint8)
    n_features = 0
    for kmer in product('ATGC', repeat=kmer_len):
        encoded = _encode_kmer(kmer)
        if lookup[encoded] != 255:
            continue
        reverse_complement = tuple(_BASE_COMPLEMENT[base] for base in reversed(kmer))
        lookup[encoded] = n_features
        lookup[_encode_kmer(reverse_complement)] = n_features
        n_features += 1

    if n_features != 136 or np.any(lookup == 255):
        raise RuntimeError(f'Invalid canonical 4-mer lookup: features={n_features}.')
    lookup.setflags(write=False)
    return lookup, n_features


def count_kmers(sequence, lookup, n_features, kmer_len=4):
    if len(sequence) < kmer_len:
        return np.zeros(n_features, dtype='<u4')

    n_windows = len(sequence) - kmer_len + 1
    if n_windows > np.iinfo(np.uint32).max:
        raise OverflowError('The sequence contains too many k-mer windows for uint32 counts.')

    sequence_bytes = np.frombuffer(sequence.encode('ascii', errors='replace'), dtype=np.uint8)
    base_codes = _ASCII_TO_BASE[sequence_bytes]
    windows = np.lib.stride_tricks.sliding_window_view(base_codes, kmer_len)
    valid_windows = windows[np.all(windows < 4, axis=1)]
    encoded = np.zeros(valid_windows.shape[0], dtype=np.uint8)
    for position in range(kmer_len):
        encoded *= 4
        encoded += valid_windows[:, position]

    feature_ids = lookup[encoded]
    counts = np.bincount(feature_ids, minlength=n_features)
    return counts.astype('<u4', copy=False)
