import os
from pathlib import Path

import numpy as np


CLUSTER_LABEL_DTYPE = np.dtype('<u4')


def kmeans_result_stem(bin_number, prefix='weight_seed_kmeans'):
    return f'{prefix}_k_{bin_number}_result'


def kmeans_label_name(bin_number, prefix='weight_seed_kmeans'):
    return kmeans_result_stem(bin_number, prefix) + '.labels.npy'


def kmeans_bins_name(bin_number, prefix='weight_seed_kmeans'):
    return kmeans_result_stem(bin_number, prefix) + '_bins'


def _expected_row_count(expected_rows):
    if isinstance(expected_rows, bool) or not isinstance(expected_rows, int) or expected_rows < 1:
        raise ValueError(f'Expected cluster-label row count must be a positive integer: {expected_rows!r}.')
    return expected_rows


def save_cluster_labels_atomic(path, labels, expected_rows):
    path = Path(path)
    expected_rows = _expected_row_count(expected_rows)
    values = np.asarray(labels)
    if values.ndim != 1 or values.shape != (expected_rows,):
        raise ValueError(
            f'Invalid cluster-label shape: file={path}, shape={values.shape}, '
            f'expected=({expected_rows},).')
    if not np.issubdtype(values.dtype, np.integer) or np.issubdtype(values.dtype, np.bool_):
        raise ValueError(f'Cluster labels must be integers: file={path}, dtype={values.dtype}.')
    if np.any(values < 0) or int(values.max()) > np.iinfo(np.uint32).max:
        raise ValueError(f'Cluster labels cannot be represented as uint32: {path}')

    path.parent.mkdir(parents=True, exist_ok=True)
    values = values.astype(CLUSTER_LABEL_DTYPE, copy=False)
    temporary_path = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    temporary_path.unlink(missing_ok=True)
    try:
        with temporary_path.open('wb') as output_handle:
            np.save(output_handle, values, allow_pickle=False)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return expected_rows


def validate_cluster_labels(path, expected_rows):
    path = Path(path)
    expected_rows = _expected_row_count(expected_rows)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty cluster-label file: {path}')
    try:
        labels = np.load(path, mmap_mode='r', allow_pickle=False)
    except (OSError, ValueError) as error:
        raise RuntimeError(f'Cannot load cluster-label file {path}: {error}') from error
    if labels.dtype != CLUSTER_LABEL_DTYPE or labels.shape != (expected_rows,):
        raise RuntimeError(
            f'Invalid cluster-label file: file={path}, dtype={labels.dtype}, '
            f'shape={labels.shape}, expected_dtype=uint32, expected_shape=({expected_rows},).')
    return labels


def _parse_legacy_label(value):
    numeric = value[5:] if value.startswith('group') else value
    if not numeric.isdecimal():
        raise ValueError(f'Unsupported legacy cluster label: {value!r}.')
    label = int(numeric)
    if label > np.iinfo(np.uint32).max:
        raise ValueError(f'Legacy cluster label is larger than uint32: {value!r}.')
    return label


def read_legacy_cluster_labels(path, expected_ids=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty legacy cluster result: {path}')
    identifiers = [] if expected_ids is None else None
    if expected_ids is not None and not isinstance(expected_ids, tuple):
        expected_ids = tuple(str(value) for value in expected_ids)
    values = [] if expected_ids is None else np.empty(len(expected_ids), dtype=CLUSTER_LABEL_DTYPE)
    row_count = 0
    seen_ids = set()
    try:
        with path.open(encoding='utf-8') as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 2 or not fields[0] or not fields[1]:
                    raise RuntimeError(
                        f'Malformed legacy cluster result: file={path}, line={line_number}, '
                        'expected two nonempty tab-separated fields.')
                contig_id, cluster_id = fields
                row_index = row_count
                if expected_ids is not None:
                    if row_index >= len(expected_ids):
                        raise RuntimeError(
                            f'Legacy cluster result has more than {len(expected_ids)} rows: {path}')
                    if contig_id != expected_ids[row_index]:
                        raise RuntimeError(
                            f'Legacy cluster identifier mismatch: file={path}, row={row_index + 1}, '
                            f'received={contig_id!r}, expected={expected_ids[row_index]!r}.')
                else:
                    if contig_id in seen_ids:
                        raise RuntimeError(
                            f'Legacy cluster result contains duplicate identifier {contig_id!r}: {path}')
                    seen_ids.add(contig_id)
                    identifiers.append(contig_id)
                try:
                    label = _parse_legacy_label(cluster_id)
                except ValueError as error:
                    raise RuntimeError(f'{error} File: {path}, line={line_number}.') from error
                if expected_ids is None:
                    values.append(label)
                else:
                    values[row_index] = label
                row_count += 1
    except UnicodeDecodeError as error:
        raise RuntimeError(f'Cannot decode legacy cluster result {path}: {error}') from error
    except OSError as error:
        raise RuntimeError(f'Cannot read legacy cluster result {path}: {error}') from error

    if row_count == 0:
        raise RuntimeError(f'Missing or empty legacy cluster result: {path}')
    if expected_ids is not None and row_count != len(expected_ids):
        raise RuntimeError(
            f'Invalid legacy cluster row count: file={path}, rows={row_count}, '
            f'expected={len(expected_ids)}.')
    return np.asarray(values, dtype=CLUSTER_LABEL_DTYPE), \
        expected_ids if expected_ids is not None else tuple(identifiers)
