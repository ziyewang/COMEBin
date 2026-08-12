from dataclasses import dataclass
import json
import os
from pathlib import Path

import numpy as np


FEATURE_SCHEMA_VERSION = 1
KMER_FEATURES = 136


@dataclass(frozen=True)
class AugmentationBundle:
    data_path: Path
    manifest: dict
    contig_ids: tuple
    contig_lengths: np.ndarray
    augmentation_intervals: np.ndarray
    kmer_counts: np.ndarray


@dataclass(frozen=True)
class BinaryFeatureBundle:
    augmentation: AugmentationBundle
    coverage_manifest: dict
    bam_ids: tuple
    coverage_mean: np.ndarray
    coverage_var: np.ndarray

    @property
    def data_path(self):
        return self.augmentation.data_path

    @property
    def feature_manifest(self):
        return self.augmentation.manifest

    @property
    def contig_ids(self):
        return self.augmentation.contig_ids

    @property
    def contig_lengths(self):
        return self.augmentation.contig_lengths

    @property
    def augmentation_intervals(self):
        return self.augmentation.augmentation_intervals

    @property
    def kmer_counts(self):
        return self.augmentation.kmer_counts


def publish_json_atomic(final_path, payload):
    final_path = Path(final_path)
    temporary_path = final_path.with_name(f'.{final_path.name}.{os.getpid()}.tmp')
    try:
        with temporary_path.open('w', encoding='utf-8') as output_handle:
            json.dump(payload, output_handle, indent=2)
            output_handle.write('\n')
        os.replace(temporary_path, final_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def publish_ids_atomic(final_path, identifiers):
    final_path = Path(final_path)
    temporary_path = final_path.with_name(f'.{final_path.name}.{os.getpid()}.tmp')
    try:
        with temporary_path.open('w', encoding='utf-8') as output_handle:
            for identifier in identifiers:
                output_handle.write(str(identifier) + '\n')
        os.replace(temporary_path, final_path)
    finally:
        temporary_path.unlink(missing_ok=True)


def read_ids(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'{label} identifier file is missing or empty: {path}')

    with path.open(encoding='utf-8') as input_handle:
        identifiers = tuple(line.rstrip('\r\n') for line in input_handle)
    if not identifiers or any(not identifier for identifier in identifiers):
        raise ValueError(f'{label} identifier file contains an empty identifier: {path}')
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(f'{label} identifier file contains duplicate identifiers: {path}')
    return identifiers


def load_json_manifest(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'{label} manifest is missing or empty: {path}')
    try:
        with path.open(encoding='utf-8') as input_handle:
            manifest = json.load(input_handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'Invalid {label} manifest {path}: {error}') from error
    if not isinstance(manifest, dict):
        raise ValueError(f'{label} manifest must contain a JSON object: {path}')
    if manifest.get('schema_version') != FEATURE_SCHEMA_VERSION:
        raise ValueError(
            f'Unsupported {label} schema in {path}: {manifest.get("schema_version")!r}.')
    if manifest.get('complete') is not True:
        raise ValueError(f'{label} manifest is not marked complete: {path}')
    return manifest


def load_array(path, dtype, shape, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f'{label} array is missing or empty: {path}')
    try:
        values = np.load(path, mmap_mode='r', allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot load {label} array {path}: {error}') from error
    expected_dtype = np.dtype(dtype)
    expected_shape = tuple(shape)
    if values.dtype != expected_dtype or values.shape != expected_shape:
        raise ValueError(
            f'Invalid {label} array: file={path}, dtype={values.dtype}, shape={values.shape}, '
            f'expected_dtype={expected_dtype}, expected_shape={expected_shape}.')
    return values


def _manifest_int(manifest, key, label, minimum=0):
    value = manifest.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'Invalid {key} in {label} manifest: {value!r}.')
    return value


def _validate_payload(manifest, key, filename, dtype, shape, label):
    payload = manifest.get(key)
    expected = {'file': filename, 'dtype': np.dtype(dtype).name, 'shape': list(shape)}
    if payload != expected:
        raise ValueError(
            f'Invalid {key} declaration in {label} manifest: received={payload!r}, expected={expected!r}.')


def load_augmentation_bundle(data_path, n_views=None):
    data_path = Path(data_path).expanduser().resolve()
    manifest = load_json_manifest(data_path / 'feature_manifest.json', 'Feature')
    contig_count = _manifest_int(manifest, 'contig_count', 'feature', minimum=1)
    manifest_views = _manifest_int(manifest, 'n_views', 'feature', minimum=1)
    if n_views is not None and manifest_views != n_views:
        raise ValueError(
            f'Feature view count mismatch: manifest={manifest_views}, requested={n_views}.')
    assembly_sequence_count = _manifest_int(
        manifest, 'assembly_sequence_count', 'feature', minimum=1)
    assembly_n50 = _manifest_int(manifest, 'assembly_n50', 'feature', minimum=1)
    contig_threshold = _manifest_int(
        manifest, 'contig_length_threshold', 'feature', minimum=0)
    if assembly_sequence_count < contig_count:
        raise ValueError(
            f'Feature manifest has fewer assembly sequences than retained contigs: '
            f'assembly={assembly_sequence_count}, retained={contig_count}.')
    if manifest.get('interval_convention') != 'zero-based-half-open':
        raise ValueError(
            f'Unsupported interval convention: {manifest.get("interval_convention")!r}.')

    identifiers = manifest.get('identifiers')
    expected_identifiers = {'file': 'contig_ids.txt', 'rows': contig_count}
    if identifiers != expected_identifiers:
        raise ValueError(
            f'Invalid identifier declaration in feature manifest: '
            f'received={identifiers!r}, expected={expected_identifiers!r}.')
    _validate_payload(
        manifest, 'contig_lengths', 'contig_lengths.npy', np.int64,
        (contig_count,), 'feature')
    _validate_payload(
        manifest, 'augmentation_intervals', 'augmentation_intervals.npy', np.int64,
        (contig_count, manifest_views, 2), 'feature')
    _validate_payload(
        manifest, 'kmer_counts', 'kmer_counts.npy', np.uint32,
        (contig_count, manifest_views, KMER_FEATURES), 'feature')

    contig_ids = read_ids(data_path / 'contig_ids.txt', 'Contig')
    if len(contig_ids) != contig_count:
        raise ValueError(
            f'Contig identifier count mismatch: file={len(contig_ids)}, manifest={contig_count}.')
    contig_lengths = load_array(
        data_path / 'contig_lengths.npy', np.int64, (contig_count,), 'contig length')
    intervals = load_array(
        data_path / 'augmentation_intervals.npy', np.int64,
        (contig_count, manifest_views, 2), 'augmentation interval')
    kmer_counts = load_array(
        data_path / 'kmer_counts.npy', np.uint32,
        (contig_count, manifest_views, KMER_FEATURES), 'k-mer count')

    if np.any(contig_lengths <= contig_threshold):
        raise ValueError(
            f'At least one retained contig is not longer than the manifest threshold '
            f'{contig_threshold}.')
    starts = intervals[:, :, 0]
    ends = intervals[:, :, 1]
    if np.any(starts[:, 0] != 0) or np.any(ends[:, 0] != contig_lengths):
        raise ValueError('View 0 intervals do not cover their complete contigs.')
    if np.any(starts < 0) or np.any(ends <= starts) or np.any(ends > contig_lengths[:, None]):
        raise ValueError('At least one augmentation interval is outside its reference contig.')

    return AugmentationBundle(
        data_path, manifest, contig_ids, contig_lengths, intervals, kmer_counts)


def load_binary_feature_bundle(data_path, n_views=None):
    augmentation = load_augmentation_bundle(data_path, n_views=n_views)
    data_path = augmentation.data_path
    manifest = load_json_manifest(data_path / 'coverage_manifest.json', 'Coverage')
    contig_count = len(augmentation.contig_ids)
    manifest_contigs = _manifest_int(manifest, 'contig_count', 'coverage', minimum=1)
    manifest_views = _manifest_int(manifest, 'n_views', 'coverage', minimum=1)
    sample_count = _manifest_int(manifest, 'sample_count', 'coverage', minimum=1)
    feature_views = int(augmentation.manifest['n_views'])
    if manifest_contigs != contig_count or manifest_views != feature_views:
        raise ValueError(
            f'Coverage dimensions do not match augmentation: coverage_contigs={manifest_contigs}, '
            f'feature_contigs={contig_count}, coverage_views={manifest_views}, '
            f'feature_views={feature_views}.')

    identifiers = manifest.get('identifiers')
    expected_identifiers = {'file': 'bam_ids.txt', 'rows': sample_count}
    if identifiers != expected_identifiers:
        raise ValueError(
            f'Invalid identifier declaration in coverage manifest: '
            f'received={identifiers!r}, expected={expected_identifiers!r}.')
    shape = (sample_count, manifest_views, contig_count)
    _validate_payload(
        manifest, 'coverage_mean', 'coverage_mean.npy', np.float64, shape, 'coverage')
    _validate_payload(
        manifest, 'coverage_var', 'coverage_var.npy', np.float64, shape, 'coverage')

    bam_ids = read_ids(data_path / 'bam_ids.txt', 'BAM')
    if len(bam_ids) != sample_count:
        raise ValueError(
            f'BAM identifier count mismatch: file={len(bam_ids)}, manifest={sample_count}.')
    coverage_mean = load_array(
        data_path / 'coverage_mean.npy', np.float64, shape, 'coverage mean')
    coverage_var = load_array(
        data_path / 'coverage_var.npy', np.float64, shape, 'coverage variance')
    return BinaryFeatureBundle(
        augmentation, manifest, bam_ids, coverage_mean, coverage_var)
