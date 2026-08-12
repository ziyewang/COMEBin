from dataclasses import dataclass
from itertools import groupby
import multiprocessing
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from feature_bundle import FEATURE_SCHEMA_VERSION, load_augmentation_bundle
from feature_bundle import publish_ids_atomic, publish_json_atomic


@dataclass(frozen=True)
class CoverageContext:
    contig_names: tuple
    contig_index: dict
    contig_lengths: np.ndarray
    starts: np.ndarray
    ends: np.ndarray


_COVERAGE_CONTEXT = None


def _build_coverage_context(out_path, n_views, contig_len):
    bundle = load_augmentation_bundle(out_path, n_views=n_views)
    if np.any(bundle.contig_lengths <= contig_len):
        raise ValueError(
            f'At least one retained contig is not longer than {contig_len}.')
    contig_index = {
        name: index for index, name in enumerate(bundle.contig_ids)
    }
    return CoverageContext(
        bundle.contig_ids, contig_index, bundle.contig_lengths,
        bundle.augmentation_intervals[:, :, 0],
        bundle.augmentation_intervals[:, :, 1])


def _update_stats(count, mean, m2, value, weight):
    if weight <= 0:
        return count, mean, m2
    new_count = count + weight
    delta = value - mean
    mean += delta * weight / new_count
    m2 += delta * delta * count * weight / new_count
    return new_count, mean, m2


def _parse_bedgraph_line(contig_name, line):
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) != 4 or fields[0] != contig_name:
        raise ValueError(
            f'Malformed bedGraph row for {contig_name}: {line.rstrip()}')
    try:
        depth_start = int(fields[1])
        depth_end = int(fields[2])
        depth = float(fields[3])
    except ValueError as error:
        raise ValueError(
            f'Malformed bedGraph values for {contig_name}: {line.rstrip()}') from error
    if (depth_start < 0 or depth_end <= depth_start or
            not np.isfinite(depth) or depth < 0):
        raise ValueError(
            f'Invalid bedGraph interval for {contig_name}: {line.rstrip()}')
    return depth_start, depth_end, depth


def _summarize_contig(contig_name, lines, starts, ends):
    counts = np.zeros(len(starts), dtype=np.int64)
    means = np.zeros(len(starts), dtype=np.float64)
    m2 = np.zeros(len(starts), dtype=np.float64)
    previous_end = 0
    contig_length = int(ends[0])

    for line in lines:
        depth_start, depth_end, depth = _parse_bedgraph_line(contig_name, line)
        if depth_start != previous_end or depth_end > contig_length:
            raise ValueError(
                f'Non-contiguous bedGraph coverage for {contig_name}: '
                f'previous_end={previous_end}, row={line.rstrip()}')
        previous_end = depth_end
        for view in range(len(starts)):
            overlap = max(
                0, min(depth_end, int(ends[view])) -
                max(depth_start, int(starts[view])))
            counts[view], means[view], m2[view] = _update_stats(
                counts[view], means[view], m2[view], depth, overlap)

    if previous_end != contig_length:
        raise ValueError(
            f'Incomplete bedGraph span for {contig_name}: '
            f'observed_end={previous_end}, expected_end={contig_length}.')
    expected = ends - starts
    if not np.array_equal(counts, expected):
        raise ValueError(
            f'Incomplete coverage for {contig_name}: '
            f'observed={counts.tolist()}, expected={expected.tolist()}.')
    return means, m2 / counts


def _validate_ignored_group(contig_name, lines):
    previous_end = 0
    for line in lines:
        depth_start, depth_end, _ = _parse_bedgraph_line(contig_name, line)
        if depth_start != previous_end:
            raise ValueError(
                f'Non-contiguous bedGraph coverage for ignored contig '
                f'{contig_name}: previous_end={previous_end}, row={line.rstrip()}')
        previous_end = depth_end


def _stream_bam_coverage(bam_path, context, sample_mean, sample_var, seen):
    command = ['bedtools', 'genomecov', '-bga', '-ibam', str(bam_path)]
    with tempfile.TemporaryFile(mode='w+t', encoding='utf-8') as stderr_handle:
        with subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=stderr_handle,
                text=True, encoding='utf-8', bufsize=1024 * 1024) as process:
            try:
                for contig_name, lines in groupby(
                        process.stdout, lambda line: line.split('\t', 1)[0]):
                    index = context.contig_index.get(contig_name)
                    if index is None:
                        _validate_ignored_group(contig_name, lines)
                        continue
                    if seen[index]:
                        raise ValueError(
                            f'BAM coverage contains repeated groups for retained '
                            f'contig {contig_name}: {bam_path}')
                    means, variances = _summarize_contig(
                        contig_name, lines, context.starts[index],
                        context.ends[index])
                    sample_mean[:, index] = means
                    sample_var[:, index] = variances
                    seen[index] = True
            except BaseException:
                process.terminate()
                process.wait()
                raise

            returncode = process.wait()
            if returncode != 0:
                stderr_handle.seek(0)
                stderr = stderr_handle.read()
                raise RuntimeError(
                    f'bedtools genomecov failed for {bam_path} with status '
                    f'{returncode}: {stderr[-4000:]}')


def _create_coverage_targets(out_path, n_samples, n_views, n_contigs):
    mean_path = out_path / f'.coverage_mean.{os.getpid()}.tmp.npy'
    var_path = out_path / f'.coverage_var.{os.getpid()}.tmp.npy'
    mean_path.unlink(missing_ok=True)
    var_path.unlink(missing_ok=True)
    shape = (n_samples, n_views, n_contigs)
    mean_array = np.lib.format.open_memmap(
        mean_path, mode='w+', dtype='<f8', shape=shape)
    var_array = np.lib.format.open_memmap(
        var_path, mode='w+', dtype='<f8', shape=shape)
    mean_array.flush()
    var_array.flush()
    del mean_array
    del var_array
    return mean_path, var_path


def _process_bam(job):
    bam_path, sample_index, mean_path, var_path = job
    context = _COVERAGE_CONTEXT
    if context is None:
        raise RuntimeError('Coverage worker context was not initialized.')
    expected_shape = (
        sample_index + 1, len(context.starts[0]), len(context.contig_names))
    mean_array = np.load(mean_path, mmap_mode='r+', allow_pickle=False)
    var_array = np.load(var_path, mmap_mode='r+', allow_pickle=False)
    if (mean_array.dtype != np.dtype('<f8') or
            var_array.dtype != np.dtype('<f8') or
            mean_array.shape != var_array.shape or
            len(mean_array.shape) != 3 or
            mean_array.shape[0] < expected_shape[0] or
            mean_array.shape[1:] != expected_shape[1:]):
        raise ValueError(
            f'Invalid temporary coverage targets: mean_shape={mean_array.shape}, '
            f'var_shape={var_array.shape}, sample_index={sample_index}.')

    seen = np.zeros(len(context.contig_names), dtype=bool)
    try:
        sample_mean = mean_array[sample_index]
        sample_var = var_array[sample_index]
        _stream_bam_coverage(
            bam_path, context, sample_mean, sample_var, seen)
        if not np.all(seen):
            missing = np.flatnonzero(~seen)
            examples = [context.contig_names[index] for index in missing[:5]]
            raise ValueError(
                f'BAM coverage is missing {len(missing)} retained contigs for '
                f'{bam_path}; examples={examples}.')
        if (not np.isfinite(sample_mean).all() or
                np.any(sample_mean < 0)):
            raise ValueError(f'Coverage mean contains invalid values for {bam_path}.')
        if (not np.isfinite(sample_var).all() or
                np.any(sample_var < 0)):
            raise ValueError(
                f'Coverage variance contains invalid values for {bam_path}.')
        mean_array.flush()
        var_array.flush()
    finally:
        del mean_array
        del var_array
    return sample_index


def _publish_coverage_manifest(out_path, n_contigs, n_views, n_samples):
    shape = [n_samples, n_views, n_contigs]
    manifest = {
        'schema_version': FEATURE_SCHEMA_VERSION,
        'complete': True,
        'contig_count': n_contigs,
        'n_views': n_views,
        'sample_count': n_samples,
        'identifiers': {'file': 'bam_ids.txt', 'rows': n_samples},
        'coverage_mean': {
            'file': 'coverage_mean.npy', 'dtype': 'float64',
            'shape': shape
        },
        'coverage_var': {
            'file': 'coverage_var.npy', 'dtype': 'float64',
            'shape': shape
        }
    }
    publish_json_atomic(out_path / 'coverage_manifest.json', manifest)


def run_gen_coverage(logger, args):
    if args.num_threads < 1:
        raise ValueError('num_threads must be greater than zero.')
    out_path = Path(args.out_augdata_path).expanduser().resolve()
    bam_dir = Path(args.bam_file_path).expanduser().resolve()
    if not bam_dir.is_dir():
        raise FileNotFoundError(f'BAM directory does not exist: {bam_dir}')
    bam_files = sorted(
        path.resolve() for path in bam_dir.glob('*.bam') if path.is_file())
    if not bam_files:
        raise FileNotFoundError(f'No BAM files were found in {bam_dir}.')

    coverage_manifest = out_path / 'coverage_manifest.json'
    coverage_manifest.unlink(missing_ok=True)
    context = _build_coverage_context(
        out_path, args.n_views, args.contig_len)
    workers = min(args.num_threads, len(bam_files))
    if workers < args.num_threads:
        logger.warning(
            'Coverage requested %d workers but only %d BAM files exist; '
            'using %d workers.', args.num_threads, len(bam_files), workers)
    logger.info(
        'Coverage inputs: contigs=%d, views=%d, BAMs=%d, workers=%d.',
        len(context.contig_names), args.n_views, len(bam_files), workers)

    mean_path, var_path = _create_coverage_targets(
        out_path, len(bam_files), args.n_views,
        len(context.contig_names))
    global _COVERAGE_CONTEXT
    _COVERAGE_CONTEXT = context
    jobs = [
        (bam_path, sample_index, mean_path, var_path)
        for sample_index, bam_path in enumerate(bam_files)
    ]
    try:
        with multiprocessing.get_context('fork').Pool(workers) as pool:
            pending = [
                pool.apply_async(_process_bam, (job,)) for job in jobs
            ]
            completed = [result.get() for result in pending]
        if sorted(completed) != list(range(len(bam_files))):
            raise RuntimeError(
                f'Coverage workers returned an incomplete sample set: '
                f'{sorted(completed)}.')

        os.replace(mean_path, out_path / 'coverage_mean.npy')
        os.replace(var_path, out_path / 'coverage_var.npy')
        publish_ids_atomic(
            out_path / 'bam_ids.txt',
            [bam_path.name for bam_path in bam_files])
        _publish_coverage_manifest(
            out_path, len(context.contig_names), args.n_views,
            len(bam_files))
        logger.info(
            'Published binary coverage features: BAMs=%d, views=%d, '
            'contigs=%d, dtype=float64.',
            len(bam_files), args.n_views, len(context.contig_names))
    except BaseException:
        mean_path.unlink(missing_ok=True)
        var_path.unlink(missing_ok=True)
        raise
    finally:
        _COVERAGE_CONTEXT = None
