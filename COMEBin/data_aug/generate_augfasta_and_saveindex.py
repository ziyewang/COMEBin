from array import array
import os
from pathlib import Path
import random
import shutil
import tempfile

from Bio import SeqIO
import numpy as np

from feature_bundle import FEATURE_SCHEMA_VERSION, publish_json_atomic
from sequence_stats import calculate_n50
from .gen_kmer import build_feature_lookup, count_kmers


def _open_uncompressed_fasta(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'Input FASTA is missing or empty: {path}')
    input_handle = path.open('rt', encoding='utf-8')
    first_character = input_handle.read(1)
    input_handle.seek(0)
    if first_character != '>':
        input_handle.close()
        raise ValueError(
            f'Canonical COMEBin input must be an uncompressed FASTA file: {path}')
    return input_handle


def _publish_raw_array(raw_path, final_path, dtype, shape, chunk_rows=8192):
    raw_path = Path(raw_path)
    final_path = Path(final_path)
    expected_bytes = int(np.prod(shape, dtype=np.int64)) * np.dtype(dtype).itemsize
    actual_bytes = raw_path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(
            f'Raw array size mismatch for {raw_path}: bytes={actual_bytes}, '
            f'expected={expected_bytes}, shape={shape}, dtype={np.dtype(dtype)}.')

    temporary_path = final_path.with_name(
        f'.{final_path.name}.{os.getpid()}.tmp.npy')
    try:
        source = np.memmap(raw_path, mode='r', dtype=dtype, shape=shape)
        target = np.lib.format.open_memmap(
            temporary_path, mode='w+', dtype=dtype, shape=shape)
        try:
            for start in range(0, shape[0], chunk_rows):
                end = min(start + chunk_rows, shape[0])
                target[start:end] = source[start:end]
            target.flush()
        finally:
            del target
            del source
        os.replace(temporary_path, final_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def _publish_feature_manifest(output_path, contig_count, n_views, contig_len,
                              assembly_sequence_count, assembly_n50):
    manifest = {
        'schema_version': FEATURE_SCHEMA_VERSION,
        'complete': True,
        'contig_count': contig_count,
        'assembly_sequence_count': assembly_sequence_count,
        'assembly_n50': assembly_n50,
        'n_views': n_views,
        'contig_length_threshold': contig_len,
        'interval_convention': 'zero-based-half-open',
        'identifiers': {'file': 'contig_ids.txt', 'rows': contig_count},
        'contig_lengths': {
            'file': 'contig_lengths.npy', 'dtype': 'int64',
            'shape': [contig_count]
        },
        'augmentation_intervals': {
            'file': 'augmentation_intervals.npy', 'dtype': 'int64',
            'shape': [contig_count, n_views, 2]
        },
        'kmer_counts': {
            'file': 'kmer_counts.npy', 'dtype': 'uint32',
            'shape': [contig_count, n_views, 136]
        }
    }
    publish_json_atomic(output_path / 'feature_manifest.json', manifest)


def run_gen_augfasta(logger, args):
    contig_file = Path(args.contig_file).expanduser().resolve()
    output_path = Path(args.out_augdata_path).expanduser().resolve()
    n_views = args.n_views
    contig_len = args.contig_len
    if n_views < 1:
        raise ValueError('n_views must be greater than zero.')
    if contig_len < 1:
        raise ValueError('contig_len must be greater than zero.')

    output_path.mkdir(parents=True, exist_ok=True)
    manifest_path = output_path / 'feature_manifest.json'
    manifest_path.unlink(missing_ok=True)
    (output_path / 'coverage_manifest.json').unlink(missing_ok=True)
    work_path = Path(tempfile.mkdtemp(prefix='.augmentation-', dir=output_path))
    ids_path = work_path / 'contig_ids.txt'
    lengths_raw_path = work_path / 'contig_lengths.raw'
    intervals_raw_path = work_path / 'augmentation_intervals.raw'
    kmer_raw_path = work_path / 'kmer_counts.raw'
    kmer_lookup, kmer_features = build_feature_lookup(4)
    assembly_lengths = array('q')
    if assembly_lengths.itemsize != np.dtype('<i8').itemsize:
        raise RuntimeError('The platform signed-long-long width is not 64 bits.')

    view_rngs = []
    for view in range(1, n_views):
        if args.seed is None:
            view_rngs.append(random.Random())
            logger.info('Generate augmentation view %d/%d without a fixed seed.',
                        view, n_views - 1)
        else:
            view_seed = args.seed + view
            view_rngs.append(random.Random(view_seed))
            logger.info('Generate augmentation view %d/%d with seed %d.',
                        view, n_views - 1, view_seed)

    contig_count = 0
    seen_ids = set()
    try:
        with _open_uncompressed_fasta(contig_file) as input_handle:
            with ids_path.open('w', encoding='utf-8') as ids_handle, \
                    lengths_raw_path.open('wb') as lengths_handle, \
                    intervals_raw_path.open('wb') as intervals_handle, \
                    kmer_raw_path.open('wb') as kmer_handle:
                for record in SeqIO.parse(input_handle, 'fasta'):
                    contig_id = record.id
                    if not contig_id:
                        raise ValueError(f'FASTA contains an empty contig identifier: {contig_file}')
                    if contig_id in seen_ids:
                        raise ValueError(f'Duplicate contig identifier: {contig_id}')
                    seen_ids.add(contig_id)

                    sequence = str(record.seq).upper()
                    sequence_length = len(sequence)
                    assembly_lengths.append(sequence_length)
                    if sequence_length <= contig_len:
                        continue

                    intervals = np.empty((n_views, 2), dtype='<i8')
                    kmer_counts = np.empty((n_views, kmer_features), dtype='<u4')
                    intervals[0] = (0, sequence_length)
                    kmer_counts[0] = count_kmers(
                        sequence, kmer_lookup, kmer_features, 4)
                    for view, rng in enumerate(view_rngs, start=1):
                        start = rng.randint(
                            0, sequence_length - (contig_len + 1))
                        simulated_length = rng.randint(
                            contig_len, sequence_length - start)
                        end = start + simulated_length
                        intervals[view] = (start, end)
                        kmer_counts[view] = count_kmers(
                            sequence[start:end], kmer_lookup, kmer_features, 4)

                    ids_handle.write(contig_id + '\n')
                    np.asarray(sequence_length, dtype='<i8').tofile(lengths_handle)
                    intervals.tofile(intervals_handle)
                    kmer_counts.tofile(kmer_handle)
                    contig_count += 1

        if not assembly_lengths:
            raise ValueError(f'No sequences were found in {contig_file}.')
        assembly_n50 = calculate_n50(np.asarray(assembly_lengths, dtype=np.int64))
        if contig_count == 0:
            raise ValueError(
                f'No contigs longer than {contig_len} were found in {contig_file}.')

        _publish_raw_array(
            lengths_raw_path, output_path / 'contig_lengths.npy',
            '<i8', (contig_count,))
        _publish_raw_array(
            intervals_raw_path, output_path / 'augmentation_intervals.npy',
            '<i8', (contig_count, n_views, 2))
        _publish_raw_array(
            kmer_raw_path, output_path / 'kmer_counts.npy',
            '<u4', (contig_count, n_views, kmer_features))
        os.replace(ids_path, output_path / 'contig_ids.txt')
        _publish_feature_manifest(
            output_path, contig_count, n_views, contig_len,
            len(assembly_lengths), assembly_n50)
        logger.info(
            'Generated binary augmentation features: assembly_contigs=%d, '
            'assembly_n50=%d, retained_contigs=%d, views=%d, kmer_features=%d.',
            len(assembly_lengths), assembly_n50, contig_count, n_views,
            kmer_features)
    finally:
        shutil.rmtree(work_path, ignore_errors=True)
