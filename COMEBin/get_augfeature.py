import logging
from pathlib import Path

import numpy as np
import torch

from feature_bundle import load_binary_feature_bundle


def _normalise_sample(values, mode, label, logger):
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f'{label} contains a non-finite or negative value.')
    values += 1e-5
    if mode == 'mean':
        offset = 0.0
        denominator = values.mean()
    elif mode == 'minmax':
        offset = values.min()
        denominator = values.max() - offset
    elif mode == 'standard':
        offset = values.mean()
        denominator = values.std()
    else:
        offset = 0.0
        denominator = values.max()

    if not np.isfinite(denominator) or denominator < 0:
        raise ValueError(
            f'{label} has an invalid normalization denominator: {denominator}.')
    if denominator == 0:
        logger.warning(
            '%s is constant; its normalized values are set to zero.', label)
        values.fill(0.0)
    else:
        values -= offset
        values /= denominator


def get_kmer_coverage(
        data_path, n_views=2, cov_meannormalize=False,
        cov_minmaxnormalize=False, cov_standardization=False, addvars=False,
        vars_sqrt=False, logger=None):
    logger = logger or logging.getLogger('COMEBin')
    data_path = Path(data_path).expanduser().resolve()
    bundle = load_binary_feature_bundle(data_path, n_views=n_views)
    n_contigs = len(bundle.contig_ids)
    n_samples = len(bundle.bam_ids)
    logger.info(
        'Feature assembly: validated binary bundle; contigs=%d, views=%d, '
        'BAMs=%d.', n_contigs, n_views, n_samples)

    kmer_values = np.array(
        bundle.kmer_counts.transpose(1, 0, 2), dtype=np.float64,
        order='C', copy=True)
    kmer_values += 1.0
    row_sums = kmer_values.sum(axis=2, keepdims=True)
    if not np.isfinite(row_sums).all() or np.any(row_sums <= 0):
        raise ValueError('Invalid k-mer row sum in kmer_counts.npy.')
    kmer_values /= row_sums

    variance_dimensions = n_samples if addvars else 0
    coverage_start = variance_dimensions
    kmer_start = variance_dimensions + n_samples
    feature_dimensions = variance_dimensions + n_samples + 136
    feature_tensor = torch.empty(
        (n_views, n_contigs, feature_dimensions), dtype=torch.float32)
    feature_array = feature_tensor.numpy()

    feature_array[:, :, kmer_start:] = kmer_values
    del kmer_values

    normalization_mode = 'max'
    if cov_meannormalize:
        normalization_mode = 'mean'
    elif cov_minmaxnormalize:
        normalization_mode = 'minmax'
    elif cov_standardization:
        normalization_mode = 'standard'

    for sample_index, sample_name in enumerate(bundle.bam_ids):
        values = np.array(
            bundle.coverage_mean[sample_index], dtype=np.float64, copy=True)
        _normalise_sample(
            values, normalization_mode,
            f'Coverage sample {sample_name}', logger)
        feature_array[:, :, coverage_start + sample_index] = values

        if addvars:
            values = np.array(
                bundle.coverage_var[sample_index],
                dtype=np.float64, copy=True)
            if not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(
                    f'Variance sample {sample_name} contains a non-finite '
                    f'or negative value.')
            if vars_sqrt:
                np.sqrt(values, out=values)
            _normalise_sample(
                values, normalization_mode,
                f'Variance sample {sample_name}', logger)
            feature_array[:, :, sample_index] = values

    if not np.isfinite(feature_array).all():
        raise ValueError(
            f'Invalid assembled feature tensor: shape={feature_array.shape}.')
    logger.info(
        'Feature assembly: binary tensor ready; shape=%s, dtype=%s, '
        'column_order=variance|coverage|kmer.',
        feature_array.shape, feature_array.dtype)
    return (
        list(feature_tensor.unbind(0)),
        np.asarray(bundle.contig_ids),
        bundle.contig_lengths,
        dict(bundle.feature_manifest))


def get_ContrastiveLearningDataset(
        data_path, n_views=2, cov_meannormalize=False,
        cov_minmaxnormalize=False, cov_standardization=False, addvars=False,
        vars_sqrt=False, logger=None):
    return get_kmer_coverage(
        data_path, n_views, cov_meannormalize, cov_minmaxnormalize,
        cov_standardization, addvars, vars_sqrt, logger)


def load_nocontrast_view0(data_path, logger):
    bundle = load_binary_feature_bundle(data_path)
    n_contigs = len(bundle.contig_ids)
    n_samples = len(bundle.bam_ids)
    coverage = np.empty((n_contigs, n_samples), dtype=np.float32)
    for sample_index, sample_name in enumerate(bundle.bam_ids):
        values = np.array(
            bundle.coverage_mean[sample_index, 0, :],
            dtype=np.float64, copy=True)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(
                f'Invalid view-0 coverage values for BAM {sample_name}.')
        values += 1e-5
        denominator = values.max()
        if not np.isfinite(denominator) or denominator <= 0:
            raise ValueError(
                f'Invalid view-0 coverage maximum for BAM '
                f'{sample_name}: {denominator}.')
        coverage[:, sample_index] = values / denominator

    kmer_work = np.array(
        bundle.kmer_counts[:, 0, :], dtype=np.float64,
        order='C', copy=True)
    kmer_work += 1.0
    row_sums = kmer_work.sum(axis=1, keepdims=True)
    if not np.isfinite(row_sums).all() or np.any(row_sums <= 0):
        raise ValueError('Invalid view-0 k-mer row sum.')
    kmer_work /= row_sums
    kmer = np.asarray(kmer_work, dtype=np.float32)
    del kmer_work

    if not np.isfinite(coverage).all() or not np.isfinite(kmer).all():
        raise ValueError('Non-finite value in assembled nocontrast features.')
    logger.info(
        'Loaded binary nocontrast features: contigs=%d, '
        'coverage_dimensions=%d, kmer_dimensions=%d, dtype=float32.',
        n_contigs, coverage.shape[1], kmer.shape[1])
    return (
        coverage, kmer, np.asarray(bundle.contig_ids),
        bundle.contig_lengths)
