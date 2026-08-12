from pathlib import Path

from cluster_labels import validate_cluster_labels


LEIDEN_RESOLUTIONS = [1, 5, 10, 30, 50, 70, 90, 110]
LEIDEN_BANDWIDTHS = [0.05, 0.1, 0.15, 0.2, 0.3]
LEIDEN_PARTGRAPH_RATIOS = [50, 100, 80]


def leiden_output_stem(max_edges, partgraph_ratio, bandwidth, resolution):
    return 'Leiden_bandwidth_' + str(bandwidth) + '_res_maxedges' + str(max_edges) + \
           'respara_' + str(resolution) + '_partgraph_ratio_' + str(partgraph_ratio)


def leiden_output_name(max_edges, partgraph_ratio, bandwidth, resolution):
    return leiden_output_stem(max_edges, partgraph_ratio, bandwidth, resolution) + '.labels.npy'


def legacy_leiden_output_name(max_edges, partgraph_ratio, bandwidth, resolution):
    return leiden_output_stem(max_edges, partgraph_ratio, bandwidth, resolution) + '.tsv'


def expected_leiden_output_names(max_edges):
    return [leiden_output_name(max_edges, ratio, bandwidth, resolution)
            for ratio in LEIDEN_PARTGRAPH_RATIOS
            for bandwidth in LEIDEN_BANDWIDTHS
            for resolution in LEIDEN_RESOLUTIONS]


def expected_legacy_leiden_output_names(max_edges):
    return [legacy_leiden_output_name(max_edges, ratio, bandwidth, resolution)
            for ratio in LEIDEN_PARTGRAPH_RATIOS
            for bandwidth in LEIDEN_BANDWIDTHS
            for resolution in LEIDEN_RESOLUTIONS]


def validate_leiden_tsv(path, expected_rows=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty Leiden result: {path}')

    row_count = 0
    try:
        with path.open(encoding='utf-8') as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 2 or not fields[0] or not fields[1]:
                    raise RuntimeError(
                        f'Malformed Leiden result: file={path}, line={line_number}, '
                        'expected two nonempty tab-separated fields.')
                row_count = line_number
    except UnicodeDecodeError as error:
        raise RuntimeError(f'Cannot decode Leiden result {path}: {error}') from error
    except OSError as error:
        raise RuntimeError(f'Cannot read Leiden result {path}: {error}') from error

    if row_count == 0:
        raise RuntimeError(f'Missing or empty Leiden result: {path}')
    if expected_rows is not None and row_count != expected_rows:
        raise RuntimeError(
            f'Invalid Leiden row count: file={path}, rows={row_count}, expected={expected_rows}.')
    return row_count


def validate_leiden_output_set(result_dir, max_edges, expected_rows=None):
    result_dir = Path(result_dir)
    if expected_rows is None:
        raise ValueError('Canonical Leiden validation requires the shared identifier count.')
    expected_names = expected_leiden_output_names(max_edges)
    missing_names = [name for name in expected_names if not (result_dir / name).is_file()]
    if missing_names:
        raise RuntimeError(
            f'Leiden result set is incomplete for max_edges={max_edges}; '
            f'missing {len(missing_names)}/{len(expected_names)} files: ' + ', '.join(missing_names))

    for name in expected_names:
        labels = validate_cluster_labels(result_dir / name, expected_rows)
        del labels
    return expected_names, expected_rows


def detect_leiden_output_set(result_dir, max_edges):
    result_dir = Path(result_dir)
    canonical_names = expected_leiden_output_names(max_edges)
    legacy_names = expected_legacy_leiden_output_names(max_edges)
    canonical_present = [name for name in canonical_names if (result_dir / name).is_file()]
    legacy_present = [name for name in legacy_names if (result_dir / name).is_file()]

    if canonical_present:
        if len(canonical_present) != len(canonical_names):
            missing = [name for name in canonical_names if name not in canonical_present]
            raise RuntimeError(
                f'Canonical Leiden result set is incomplete for max_edges={max_edges}; '
                f'missing {len(missing)}/{len(canonical_names)} files: ' + ', '.join(missing))
        return 'binary-labels', canonical_names
    if legacy_present:
        if len(legacy_present) != len(legacy_names):
            missing = [name for name in legacy_names if name not in legacy_present]
            raise RuntimeError(
                f'Legacy Leiden result set is incomplete for max_edges={max_edges}; '
                f'missing {len(missing)}/{len(legacy_names)} files: ' + ', '.join(missing))
        empty = [name for name in legacy_names if (result_dir / name).stat().st_size == 0]
        if empty:
            raise RuntimeError('Legacy Leiden result set contains empty files: ' + ', '.join(empty))
        return 'legacy-tsv', legacy_names
    raise RuntimeError(
        f'No complete Leiden result set was found for max_edges={max_edges}; '
        f'expected {len(canonical_names)} binary label vectors or {len(legacy_names)} legacy TSV files.')
