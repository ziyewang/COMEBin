import gzip
import os
from collections import OrderedDict
from pathlib import Path
import shutil
import tempfile

import numpy as np

from cluster_labels import CLUSTER_LABEL_DTYPE


def _open_fasta(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'FASTA file is missing or empty: {path}')
    if path.suffix == '.gz':
        return gzip.open(path, mode='rt', encoding='utf-8')
    return path.open(encoding='utf-8')


def _iter_fasta_records(path):
    with _open_fasta(path) as input_handle:
        record_id = None
        record_lines = []
        sequence_length = 0
        for line_number, line in enumerate(input_handle, 1):
            if line.startswith('>'):
                if record_id is not None:
                    yield record_id, ''.join(record_lines), sequence_length
                header = line[1:].strip()
                if not header:
                    raise ValueError(f'Empty FASTA identifier in {path} at line {line_number}.')
                record_id = header.split(None, 1)[0]
                record_lines = [line]
                sequence_length = 0
            else:
                if record_id is None:
                    raise ValueError(f'FASTA sequence data appears before the first header in {path}.')
                record_lines.append(line)
                sequence_length += len(''.join(line.split()))
        if record_id is None:
            raise ValueError(f'No FASTA records were found in {path}.')
        yield record_id, ''.join(record_lines), sequence_length


class BoundedBinHandles:
    def __init__(self, output_dir, max_open_files=128):
        if isinstance(max_open_files, bool) or not isinstance(max_open_files, int) or max_open_files < 1:
            raise ValueError('max_open_files must be a positive integer.')
        self.output_dir = Path(output_dir)
        self.max_open_files = max_open_files
        self.handles = OrderedDict()

    def write(self, label, record):
        label = int(label)
        handle = self.handles.pop(label, None)
        if handle is None:
            if len(self.handles) >= self.max_open_files:
                _, oldest_handle = self.handles.popitem(last=False)
                oldest_handle.close()
            handle = (self.output_dir / f'{label}.fa').open('a', encoding='utf-8')
        self.handles[label] = handle
        handle.write(record)
        if record and not record.endswith('\n'):
            handle.write('\n')

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def _validate_aligned_values(contig_ids, contig_lengths, labels):
    contig_ids = tuple(str(value) for value in contig_ids)
    contig_lengths = np.asarray(contig_lengths)
    labels = np.asarray(labels)
    if not contig_ids or any(not value for value in contig_ids) or len(set(contig_ids)) != len(contig_ids):
        raise ValueError('Bin writing requires nonempty, unique contig identifiers.')
    if contig_lengths.dtype != np.dtype('<i8') or contig_lengths.shape != (len(contig_ids),):
        raise ValueError(
            f'Invalid aligned contig lengths: dtype={contig_lengths.dtype}, '
            f'shape={contig_lengths.shape}, expected=({len(contig_ids)},).')
    if np.any(contig_lengths <= 0):
        raise ValueError('Aligned contig lengths contain a non-positive value.')
    if labels.dtype != CLUSTER_LABEL_DTYPE or labels.shape != (len(contig_ids),):
        raise ValueError(
            f'Invalid aligned cluster labels: dtype={labels.dtype}, shape={labels.shape}, '
            f'expected_dtype=uint32, expected_shape=({len(contig_ids)},).')
    return contig_ids, contig_lengths, labels


def write_bins_from_labels(fasta_path, contig_ids, contig_lengths, labels, output_dir,
                           keep_labels=None, minimum_length=None, max_open_files=128,
                           replace_existing=False):
    contig_ids, contig_lengths, labels = _validate_aligned_values(
        contig_ids, contig_lengths, labels)
    if minimum_length is not None and (
            isinstance(minimum_length, bool) or not isinstance(minimum_length, int) or minimum_length < 1):
        raise ValueError('minimum_length must be a positive integer when specified.')
    keep_labels = None if keep_labels is None else {int(value) for value in keep_labels}
    positions = {contig_id: index for index, contig_id in enumerate(contig_ids)}
    seen = np.zeros(len(contig_ids), dtype=bool)
    next_index = 0
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(tempfile.mkdtemp(prefix=f'.{output_dir.name}.', dir=output_dir.parent))
    handles = BoundedBinHandles(temporary_dir, max_open_files=max_open_files)
    try:
        for contig_id, record, sequence_length in _iter_fasta_records(fasta_path):
            index = positions.get(contig_id)
            if index is None:
                if minimum_length is not None and sequence_length >= minimum_length:
                    raise ValueError(
                        f'Unexpected retained FASTA identifier {contig_id!r}; it is absent from the cluster axis.')
                continue
            if seen[index]:
                raise ValueError(f'Duplicate retained FASTA identifier: {contig_id!r}.')
            if index != next_index:
                expected = contig_ids[next_index] if next_index < len(contig_ids) else '<end>'
                raise ValueError(
                    f'Retained FASTA identifiers are out of order: received={contig_id!r}, '
                    f'expected={expected!r}.')
            if sequence_length != int(contig_lengths[index]):
                raise ValueError(
                    f'FASTA length mismatch for {contig_id!r}: received={sequence_length}, '
                    f'expected={int(contig_lengths[index])}.')
            seen[index] = True
            next_index += 1
            label = int(labels[index])
            if keep_labels is None or label in keep_labels:
                handles.write(label, record)
        handles.close()
        if next_index != len(contig_ids):
            raise ValueError(
                f'FASTA is missing {len(contig_ids) - next_index} retained identifiers; '
                f'first missing={contig_ids[next_index]!r}.')
        if output_dir.exists():
            if not replace_existing:
                raise FileExistsError(f'Bin output directory already exists: {output_dir}')
            if not output_dir.is_dir():
                raise FileExistsError(f'Bin output path is not a directory: {output_dir}')
            shutil.rmtree(output_dir)
        os.replace(temporary_dir, output_dir)
    except BaseException:
        handles.close()
        shutil.rmtree(temporary_dir, ignore_errors=True)
        raise
    return next_index


def _grouped_row_indexes(labels):
    unique_labels, first_indexes, inverse = np.unique(
        labels, return_index=True, return_inverse=True)
    rank_by_unique = np.empty(len(unique_labels), dtype=np.int64)
    rank_by_unique[np.argsort(first_indexes)] = np.arange(len(unique_labels), dtype=np.int64)
    row_ranks = rank_by_unique[inverse]
    return np.argsort(row_ranks, kind='stable')


def save_assignments_atomic(path, contig_ids, labels, keep_labels=None, label_names=None,
                            group_by_label=False):
    path = Path(path)
    contig_ids = tuple(str(value) for value in contig_ids)
    labels = np.asarray(labels)
    if labels.dtype != CLUSTER_LABEL_DTYPE or labels.shape != (len(contig_ids),):
        raise ValueError('Assignment identifiers and uint32 labels are not aligned.')
    keep_labels = None if keep_labels is None else {int(value) for value in keep_labels}
    row_indexes = _grouped_row_indexes(labels) if group_by_label else range(len(contig_ids))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    temporary_path.unlink(missing_ok=True)
    try:
        with temporary_path.open('w', encoding='utf-8') as output_handle:
            for row_index in row_indexes:
                label = int(labels[row_index])
                if keep_labels is not None and label not in keep_labels:
                    continue
                cluster_name = f'group{label}' if label_names is None else str(label_names[label])
                output_handle.write(f'{contig_ids[row_index]}\t{cluster_name}\n')
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)
