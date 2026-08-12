from dataclasses import dataclass
import os
from pathlib import Path

import numpy as np

from feature_bundle import load_array, load_json_manifest, publish_ids_atomic
from feature_bundle import publish_json_atomic, read_ids


EMBEDDING_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class EmbeddingAxis:
    output_path: Path
    manifest: dict
    contig_ids: tuple
    contig_lengths: np.ndarray
    embedding_shape: tuple
    coverage_shape: tuple


@dataclass(frozen=True)
class EmbeddingBundle:
    output_path: Path
    manifest: dict
    contig_ids: tuple
    embeddings: np.ndarray
    coverage_embeddings: np.ndarray
    contig_lengths: np.ndarray


def _positive_int(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'Invalid {label} in embedding manifest: {value!r}.')
    return value


def _payload_shape(manifest, key, filename, dtype, rows):
    payload = manifest.get(key)
    if not isinstance(payload, dict):
        raise ValueError(f'Missing {key} declaration in embedding manifest.')
    shape = payload.get('shape')
    if (payload.get('file') != filename or payload.get('dtype') != np.dtype(dtype).name or
            not isinstance(shape, list) or len(shape) != 2 or
            shape[0] != rows or isinstance(shape[1], bool) or
            not isinstance(shape[1], int) or shape[1] < 1 or
            set(payload) != {'file', 'dtype', 'shape'}):
        raise ValueError(f'Invalid {key} declaration in embedding manifest: {payload!r}.')
    return tuple(shape)


def invalidate_embedding_manifest(output_path):
    output_path = Path(output_path)
    (output_path / 'embedding_manifest.json').unlink(missing_ok=True)


def publish_embedding_bundle(output_path, embedding_temp, coverage_temp,
                             contig_ids, contig_lengths, out_dim, coverage_dim):
    output_path = Path(output_path)
    embedding_temp = Path(embedding_temp)
    coverage_temp = Path(coverage_temp)
    n_contigs = len(contig_ids)
    out_dim = _positive_int(out_dim, 'embedding dimension')
    coverage_dim = _positive_int(coverage_dim, 'coverage embedding dimension')
    if n_contigs < 1:
        raise ValueError('Cannot publish an empty embedding bundle.')
    if len(set(contig_ids)) != n_contigs or any(not str(identifier) for identifier in contig_ids):
        raise ValueError('Embedding identifiers must be nonempty and unique.')

    lengths = np.asarray(contig_lengths)
    if lengths.dtype != np.dtype('<i8') or lengths.shape != (n_contigs,) or np.any(lengths <= 0):
        raise ValueError(
            f'Invalid aligned embedding lengths: dtype={lengths.dtype}, shape={lengths.shape}.')
    embedding_values = load_array(
        embedding_temp, np.float32, (n_contigs, out_dim), 'temporary embedding')
    coverage_values = load_array(
        coverage_temp, np.float32, (n_contigs, coverage_dim),
        'temporary coverage embedding')
    if not np.isfinite(embedding_values).all():
        raise ValueError(f'Temporary embeddings contain a non-finite value: {embedding_temp}')
    if not np.isfinite(coverage_values).all():
        raise ValueError(
            f'Temporary coverage embeddings contain a non-finite value: {coverage_temp}')
    del embedding_values, coverage_values

    invalidate_embedding_manifest(output_path)
    length_temp = output_path / f'.embedding_lengths.{os.getpid()}.tmp.npy'
    length_temp.unlink(missing_ok=True)
    try:
        length_target = np.lib.format.open_memmap(
            length_temp, mode='w+', dtype='<i8', shape=(n_contigs,))
        try:
            length_target[:] = lengths
            length_target.flush()
        finally:
            del length_target

        os.replace(embedding_temp, output_path / 'embeddings.npy')
        os.replace(coverage_temp, output_path / 'covembeddings.npy')
        os.replace(length_temp, output_path / 'embedding_lengths.npy')
        publish_ids_atomic(output_path / 'embedding_ids.txt', contig_ids)
        manifest = {
            'schema_version': EMBEDDING_SCHEMA_VERSION,
            'complete': True,
            'contig_count': n_contigs,
            'identifiers': {'file': 'embedding_ids.txt', 'rows': n_contigs},
            'embeddings': {
                'file': 'embeddings.npy', 'dtype': 'float32',
                'shape': [n_contigs, out_dim]
            },
            'coverage_embeddings': {
                'file': 'covembeddings.npy', 'dtype': 'float32',
                'shape': [n_contigs, coverage_dim]
            },
            'lengths': {
                'file': 'embedding_lengths.npy', 'dtype': 'int64',
                'shape': [n_contigs]
            }
        }
        publish_json_atomic(output_path / 'embedding_manifest.json', manifest)
    finally:
        embedding_temp.unlink(missing_ok=True)
        coverage_temp.unlink(missing_ok=True)
        length_temp.unlink(missing_ok=True)


def load_embedding_axis(path):
    path = Path(path).expanduser().resolve()
    output_path = path.parent if path.name in ('embeddings.npy', 'embedding_manifest.json') else path
    manifest = load_json_manifest(
        output_path / 'embedding_manifest.json', 'Embedding')
    if manifest.get('schema_version') != EMBEDDING_SCHEMA_VERSION:
        raise ValueError(
            f'Unsupported embedding schema: {manifest.get("schema_version")!r}.')
    n_contigs = _positive_int(manifest.get('contig_count'), 'contig count')
    identifiers = manifest.get('identifiers')
    expected_identifiers = {'file': 'embedding_ids.txt', 'rows': n_contigs}
    if identifiers != expected_identifiers:
        raise ValueError(
            f'Invalid identifier declaration in embedding manifest: {identifiers!r}.')
    embedding_shape = _payload_shape(
        manifest, 'embeddings', 'embeddings.npy', np.float32, n_contigs)
    coverage_shape = _payload_shape(
        manifest, 'coverage_embeddings', 'covembeddings.npy', np.float32,
        n_contigs)
    expected_lengths = {
        'file': 'embedding_lengths.npy', 'dtype': 'int64',
        'shape': [n_contigs]
    }
    if manifest.get('lengths') != expected_lengths:
        raise ValueError(
            f'Invalid lengths declaration in embedding manifest: '
            f'{manifest.get("lengths")!r}.')

    contig_ids = read_ids(output_path / 'embedding_ids.txt', 'Embedding')
    if len(contig_ids) != n_contigs:
        raise ValueError(
            f'Embedding identifier count mismatch: file={len(contig_ids)}, '
            f'manifest={n_contigs}.')
    contig_lengths = load_array(
        output_path / 'embedding_lengths.npy', np.int64, (n_contigs,),
        'embedding length')
    if np.any(contig_lengths <= 0):
        raise ValueError('Embedding length array contains a non-positive value.')
    return EmbeddingAxis(
        output_path, manifest, contig_ids, contig_lengths,
        embedding_shape, coverage_shape)


def load_embedding_bundle(embedding_file):
    embedding_file = Path(embedding_file).expanduser().resolve()
    axis = load_embedding_axis(embedding_file)
    output_path = axis.output_path
    if embedding_file.name != axis.manifest['embeddings']['file']:
        raise ValueError(
            f'Canonical embedding input must be {axis.manifest["embeddings"]["file"]}: '
            f'{embedding_file}')
    embeddings = load_array(
        embedding_file, np.float32, axis.embedding_shape, 'embedding')
    coverage_embeddings = load_array(
        output_path / 'covembeddings.npy', np.float32, axis.coverage_shape,
        'coverage embedding')
    if not np.isfinite(embeddings).all():
        raise ValueError(f'Embedding array contains a non-finite value: {embedding_file}')
    if not np.isfinite(coverage_embeddings).all():
        raise ValueError(
            f'Coverage embedding array contains a non-finite value: '
            f'{output_path / "covembeddings.npy"}')
    return EmbeddingBundle(
        output_path, axis.manifest, axis.contig_ids, embeddings,
        coverage_embeddings, axis.contig_lengths)
