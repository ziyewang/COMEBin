import os
from collections import defaultdict
from pathlib import Path

import numpy as np

from bin_writer import save_assignments_atomic, write_bins_from_labels
from cluster_labels import kmeans_bins_name, kmeans_label_name, read_legacy_cluster_labels
from cluster_labels import validate_cluster_labels
from embedding_bundle import load_embedding_axis
from leiden_grid import detect_leiden_output_set, validate_leiden_output_set
from scripts.unitem_markers import Markers
from scripts.unitem_profile import Profile
from utils import get_length


MINIMUM_FINAL_BIN_SIZE = 200000


def get_binstats(bin_contig_names, markers):
    _, comp, cont = markers.bin_quality(bin_contig_names)
    return comp, cont


def get_bin_quality(orig_bins, methods_sorted, markers):
    """Legacy in-memory scorer retained as an equivalence-test reference."""
    bin_quality_dict = {}
    best_method = None
    best_key = None
    for method_id in methods_sorted:
        stats = _empty_stats()
        for contig_names in orig_bins[method_id].values():
            comp, cont = get_binstats(contig_names, markers)
            _record_quality(stats, comp, cont)
        _finish_stats(stats)
        bin_quality_dict[method_id] = stats
        comparison_key = (stats['sum'], stats['sum_cont5'])
        if best_key is None or comparison_key > best_key:
            best_method = method_id
            best_key = comparison_key
    if best_method is None:
        raise RuntimeError('No clustering methods were available for quality estimation.')
    return bin_quality_dict, best_method


def _empty_stats():
    return {'num_5010': 0, 'num_7010': 0, 'num_9010': 0,
            'num_505': 0, 'num_705': 0, 'num_905': 0}


def _record_quality(stats, comp, cont):
    stats['num_5010'] += int(comp > 50 and cont < 10)
    stats['num_7010'] += int(comp > 70 and cont < 10)
    stats['num_9010'] += int(comp > 90 and cont < 10)
    stats['num_505'] += int(comp > 50 and cont < 5)
    stats['num_705'] += int(comp > 70 and cont < 5)
    stats['num_905'] += int(comp > 90 and cont < 5)


def _finish_stats(stats):
    stats['sum'] = sum(stats.values())
    stats['sum_cont5'] = stats['num_505'] + stats['num_705'] + stats['num_905']


def marker_contig_indexes(contig_ids, markers):
    marker_ids = set(markers.bac_markers_on_contigs)
    marker_ids.update(markers.ar_markers_on_contigs)
    indexes = []
    names = []
    found = set()
    for index, contig_id in enumerate(contig_ids):
        if contig_id in marker_ids:
            indexes.append(index)
            names.append(contig_id)
            found.add(contig_id)
    missing = marker_ids.difference(found)
    if missing:
        raise RuntimeError(
            f'{len(missing)} marker-bearing contigs are absent from the clustering axis; '
            f'examples={sorted(missing)[:5]}.')
    return np.asarray(indexes, dtype=np.int64), tuple(names)


def score_marker_labels(labels, marker_indexes, marker_names, markers):
    marker_bins = defaultdict(list)
    for contig_name, label in zip(marker_names, labels[marker_indexes]):
        marker_bins[int(label)].append(contig_name)
    stats = _empty_stats()
    qualities = {}
    for label, contig_names in marker_bins.items():
        comp, cont = get_binstats(contig_names, markers)
        qualities[label] = (comp, cont)
        _record_quality(stats, comp, cont)
    _finish_stats(stats)
    return stats, qualities


def write_estimated_bin_quality(bin_quality_dict, output_file):
    with Path(output_file).open('w', encoding='utf-8') as output_handle:
        output_handle.write(
            'Binning_method\tnum_5010\tnum_7010\tnum_9010\tnum_505\t'
            'num_705\tnum_905\tsum\tsum_cont5\n')
        for method_id, stats in bin_quality_dict.items():
            output_handle.write(
                f'{method_id}\t{stats["num_5010"]}\t{stats["num_7010"]}\t'
                f'{stats["num_9010"]}\t{stats["num_505"]}\t{stats["num_705"]}\t'
                f'{stats["num_905"]}\t{stats["sum"]}\t{stats["sum_cont5"]}\n')


def _canonical_result_axis(args):
    axis_source = args.emb_file if args.emb_file and Path(args.emb_file).suffix == '.npy' else args.output_path
    axis = load_embedding_axis(axis_source)
    retained = np.asarray(axis.contig_lengths) >= args.contig_len
    if retained.all():
        return axis.contig_ids, axis.contig_lengths
    contig_ids = tuple(contig_id for contig_id, keep in zip(axis.contig_ids, retained) if keep)
    return contig_ids, np.asarray(axis.contig_lengths[retained], dtype=np.int64)


def _legacy_result_axis(first_result, contig_file):
    _, contig_ids = read_legacy_cluster_labels(first_result)
    length_lookup = get_length(contig_file)
    missing = [contig_id for contig_id in contig_ids if contig_id not in length_lookup]
    if missing:
        raise RuntimeError(
            f'{len(missing)} legacy result identifiers are absent from {contig_file}; '
            f'examples={missing[:5]}.')
    return contig_ids, np.asarray([length_lookup[contig_id] for contig_id in contig_ids], dtype=np.int64)


def load_result_axis(args, result_format, result_dir, result_names):
    manifest = Path(args.output_path) / 'embedding_manifest.json'
    if manifest.is_file() or (args.emb_file and Path(args.emb_file).suffix == '.npy'):
        return _canonical_result_axis(args)
    if result_format == 'legacy-tsv':
        return _legacy_result_axis(Path(result_dir) / result_names[0], args.contig_file)
    raise RuntimeError(
        'Canonical Leiden labels require embedding_ids.txt, embedding_lengths.npy, '
        'and embedding_manifest.json in the COMEBin result directory.')


def _load_result_labels(result_format, path, contig_ids):
    if result_format == 'binary-labels':
        return validate_cluster_labels(path, len(contig_ids))
    labels, _ = read_legacy_cluster_labels(path, expected_ids=contig_ids)
    return labels


def select_best_result(result_dir, result_names, result_format, contig_ids, markers):
    marker_indexes, marker_names = marker_contig_indexes(contig_ids, markers)
    all_stats = {}
    best_name = None
    best_key = None
    best_qualities = None
    for result_name in sorted(result_names):
        labels = _load_result_labels(result_format, Path(result_dir) / result_name, contig_ids)
        stats, qualities = score_marker_labels(labels, marker_indexes, marker_names, markers)
        all_stats[result_name] = stats
        comparison_key = (stats['sum'], stats['sum_cont5'])
        if best_key is None or comparison_key > best_key:
            best_name = result_name
            best_key = comparison_key
            best_qualities = qualities
        del labels
    if best_name is None:
        raise RuntimeError('No valid Leiden result was available for final selection.')
    return all_stats, best_name, best_qualities


def _ordered_labels(labels):
    unique_labels, first_indexes = np.unique(labels, return_index=True)
    return [int(value) for value in unique_labels[np.argsort(first_indexes)]]


def save_high_quality_assignments(result_dir, best_name, contig_ids, labels, qualities):
    label_order = _ordered_labels(labels)
    thresholds = [('5010_res.txt', 50, 10), ('5005_res.txt', 50, 5)]
    for suffix, min_completeness, max_contamination in thresholds:
        selected = [label for label in label_order
                    if label in qualities and qualities[label][0] > min_completeness
                    and qualities[label][1] < max_contamination]
        label_names = {label: index for index, label in enumerate(selected)}
        save_assignments_atomic(
            Path(result_dir) / f'{best_name}{suffix}', contig_ids, labels,
            keep_labels=selected, label_names=label_names, group_by_label=True)


def retained_final_labels(labels, contig_lengths, minimum_size=MINIMUM_FINAL_BIN_SIZE):
    unique_labels, inverse = np.unique(labels, return_inverse=True)
    bin_sizes = np.zeros(len(unique_labels), dtype=np.int64)
    np.add.at(bin_sizes, inverse, contig_lengths)
    return {int(value) for value in unique_labels[bin_sizes >= minimum_size]}


def publish_final_result(logger, args, contig_ids, contig_lengths, labels, keep_labels):
    final_bins = Path(args.output_path) / 'comebin_res_bins'
    final_result = Path(args.output_path) / 'comebin_res.tsv'
    logger.info('Stream winning Leiden bins from %s.', args.contig_file)
    write_bins_from_labels(
        args.contig_file, contig_ids, contig_lengths, labels, final_bins,
        keep_labels=keep_labels, minimum_length=args.contig_len, replace_existing=True)
    save_assignments_atomic(
        final_result, contig_ids, labels, keep_labels=keep_labels, group_by_label=True)


def _profile_marker_tables(logger, args, seed_num, num_threads, result_format, res_name):
    cluster_res_path = Path(args.output_path) / 'cluster_res'
    if res_name is None:
        if result_format == 'binary-labels':
            res_name = kmeans_label_name(seed_num)
            bins_path = cluster_res_path / kmeans_bins_name(seed_num)
        else:
            res_name = f'weight_seed_kmeans_k_{seed_num}_result.tsv'
            bins_path = cluster_res_path / f'{res_name}_bins'
    elif res_name.endswith('.labels.npy'):
        bins_path = cluster_res_path / (res_name[:-len('.labels.npy')] + '_bins')
    else:
        bins_path = cluster_res_path / f'{res_name}_bins'
    if not bins_path.is_dir():
        raise RuntimeError(f'KMeans bin directory is missing: {bins_path}')

    output_dir = cluster_res_path / 'unitem_profile'
    if not output_dir.exists():
        logger.info("Run unitem profile:\t%s", seed_num)
        output_dir.mkdir(parents=True)
        Profile(num_threads).run({res_name: (str(bins_path), 'fa')}, str(output_dir))
    method_dir = output_dir / 'binning_methods' / res_name
    return (method_dir / 'checkm_bac/marker_gene_table.tsv',
            method_dir / 'checkm_ar/marker_gene_table.tsv')


def run_get_final_result(logger, args, seed_num: int, num_threads: int = 40,
                         res_name=None, ignore_kmeans_res: bool = True):
    logger.info("Seed_num:\t%s", seed_num)
    result_dir = Path(args.output_path) / 'cluster_res'
    result_format, result_names = detect_leiden_output_set(result_dir, args.max_edges)
    contig_ids, contig_lengths = load_result_axis(
        args, result_format, result_dir, result_names)
    if result_format == 'binary-labels':
        result_names, result_rows = validate_leiden_output_set(
            result_dir, args.max_edges, expected_rows=len(contig_ids))
    else:
        result_rows = len(contig_ids)
    logger.info(
        "Validated complete Leiden result grid before final selection: format=%s, "
        "max_edges=%d, files=%d, rows_per_file=%d.",
        result_format, args.max_edges, len(result_names), result_rows)

    if args.bac_mg_table and args.ar_mg_table:
        bac_mg_table = Path(args.bac_mg_table)
        ar_mg_table = Path(args.ar_mg_table)
    else:
        bac_mg_table, ar_mg_table = _profile_marker_tables(
            logger, args, seed_num, num_threads, result_format, res_name)

    markers = Markers()
    markers.marker_gene_tables(str(bac_mg_table), str(ar_mg_table))
    bin_quality_dict, best_method, best_qualities = select_best_result(
        result_dir, result_names, result_format, contig_ids, markers)
    write_estimated_bin_quality(bin_quality_dict, result_dir / 'estimate_res.txt')

    best_labels = _load_result_labels(result_format, result_dir / best_method, contig_ids)
    save_high_quality_assignments(
        result_dir, best_method, contig_ids, best_labels, best_qualities)
    keep_labels = retained_final_labels(best_labels, contig_lengths)
    logger.info('Final result:\t%s', result_dir / best_method)
    publish_final_result(
        logger, args, contig_ids, contig_lengths, best_labels, keep_labels)
    del best_labels
