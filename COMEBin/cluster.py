import hnswlib
import leidenalg
import numpy as np
import gc
import time
import os
import scipy.sparse as sp
import logging
import multiprocessing
import secrets
from pathlib import Path

from bin_writer import write_bins_from_labels
from cluster_labels import kmeans_bins_name, kmeans_label_name, save_cluster_labels_atomic
from cluster_labels import validate_cluster_labels
from embedding_bundle import load_embedding_bundle
from igraph import Graph
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import normalize
from sklearn.utils import check_random_state
from sklearn.utils.extmath import row_norms, stable_cumsum

from sequence_stats import calculate_n50
from leiden_grid import LEIDEN_BANDWIDTHS, LEIDEN_PARTGRAPH_RATIOS, LEIDEN_RESOLUTIONS
from leiden_grid import leiden_output_name, validate_leiden_output_set
from typing import List, Optional


logger = logging.getLogger('COMEBin')

logger.setLevel(logging.INFO)

# logging
formatter = logging.Formatter('%(asctime)s - %(message)s')

console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)

logger.addHandler(console_hdr)


_LEIDEN_CONTEXT = None


def fit_hnsw_index(logger, features, num_threads: int, ef: int = 100, M: int = 16, space: str = 'l2',
                   save_index_file: bool = False, random_seed: Optional[int] = None) -> hnswlib.Index:
    """
    Fit an HNSW index with the given features using the HNSWlib library; Convenience function to create HNSW graph.

    :param logger: The logger object for logging messages.
    :param features: A list of lists containing the embeddings.
    :param num_threads: The number of HNSW construction threads for an unseeded run.
    :param ef: The ef parameter to tune the HNSW algorithm (default: 100).
    :param M: The M parameter to tune the HNSW algorithm (default: 16).
    :param space: The space in which the index operates (default: 'l2').
    :param save_index_file: The path to save the HNSW index file (optional).

    :return: The HNSW index created using the given features.

    This function fits an HNSW index to the provided features, allowing efficient similarity search in high-dimensional spaces.
    """

    time_start = time.time()
    if num_threads < 1:
        raise ValueError('num_threads must be greater than zero.')
    num_elements = len(features)
    labels_index = np.arange(num_elements)
    EMBEDDING_SIZE = len(features[0])

    # Declaring index
    # possible space options are l2, cosine or ip
    p = hnswlib.Index(space=space, dim=EMBEDDING_SIZE)

    # Initing index - the maximum number of elements should be known
    hnsw_seed = random_seed if random_seed is not None else secrets.randbelow(2**31 - 1) + 1
    p.init_index(max_elements=num_elements, ef_construction=ef, M=M, random_seed=hnsw_seed)

    # Element insertion
    construction_threads = 1 if random_seed is not None else num_threads
    if random_seed is not None:
        logger.info("Build deterministic HNSW index with seed=%d and one construction thread.", random_seed)
    else:
        logger.info("Build unseeded HNSW index with %d construction threads.", construction_threads)
    p.add_items(features, labels_index, num_threads=construction_threads)

    # Controlling the recall by setting ef
    # ef should always be > k
    p.set_ef(ef)

    # If you want to save the graph to a file
    if save_index_file:
        p.save_index(save_index_file)
    time_end = time.time()
    logger.info('Time cost:\t' +str(time_end - time_start) + "s")
    return p


def _copy_seed_centers(X, seed_indexes):
    if sp.issparse(X):
        return X[seed_indexes].toarray()
    return np.array(X[seed_indexes], dtype=X.dtype, copy=True)


def prepare_fixed_seed_state(X, seed_indexes):
    seed_indexes = np.asarray(seed_indexes, dtype=np.int64)
    x_squared_norms = row_norms(X, squared=True)
    fixed_centers = _copy_seed_centers(X, seed_indexes)
    closest_dist_sq = euclidean_distances(
        fixed_centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms, squared=True)
    for center in fixed_centers[1:]:
        closest_dist_sq = np.minimum(
            closest_dist_sq,
            euclidean_distances(center[np.newaxis], X, Y_norm_squared=x_squared_norms, squared=True))
    return x_squared_norms, fixed_centers, closest_dist_sq


def complete_random_seed_init(X, n_clusters, random_state, seed_indexes, x_squared_norms,
                              fixed_centers, base_closest_dist_sq, n_local_trials=None):
    random_state = check_random_state(random_state)
    n_features = X.shape[1]
    centers = np.empty((n_clusters, n_features), dtype=X.dtype)
    centers[:len(fixed_centers)] = fixed_centers
    closest_dist_sq = np.array(base_closest_dist_sq, copy=True)
    current_pot = closest_dist_sq.sum()

    if n_local_trials is None:
        n_local_trials = 2 + int(np.log(n_clusters))

    for c in range(len(seed_indexes), n_clusters):
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq), rand_vals)
        np.clip(candidate_ids, None, closest_dist_sq.size - 1, out=candidate_ids)
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in range(n_local_trials):
            new_dist_sq = np.minimum(closest_dist_sq, distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()
            if best_candidate is None or new_pot < best_pot:
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers


def make_cached_partial_seed_init(seed_indexes):
    seed_indexes = np.asarray(seed_indexes, dtype=np.int64).copy()
    cached_state = None

    def cached_partial_seed_init(X, n_clusters, random_state):
        nonlocal cached_state
        if cached_state is None:
            cached_state = prepare_fixed_seed_state(X, seed_indexes)
        x_squared_norms, fixed_centers, base_closest_dist_sq = cached_state
        return complete_random_seed_init(
            X, n_clusters, random_state, seed_indexes, x_squared_norms,
            fixed_centers, base_closest_dist_sq)

    return cached_partial_seed_init


def fit_seed_kmeans(logger, X_mat, length_weight, bin_number, seed_indexes, random_seed):
    if len(seed_indexes) == bin_number:
        kmeans_init = _copy_seed_centers(X_mat, seed_indexes)
        kmeans_n_init = 1
        init_mode = 'full'
    else:
        kmeans_init = make_cached_partial_seed_init(seed_indexes)
        kmeans_n_init = 10
        init_mode = 'partial'

    logger.info(
        "Weighted seed-KMeans initialization: mode=%s, retained_seeds=%d, clusters=%d, n_init=%d.",
        init_mode, len(seed_indexes), bin_number, kmeans_n_init)
    fit_start = time.perf_counter()
    km = KMeans(n_clusters=bin_number, random_state=random_seed, algorithm="lloyd",
                init=kmeans_init, n_init=kmeans_n_init)
    km.fit(X_mat, sample_weight=length_weight)
    logger.info("Weighted seed-KMeans fit completed: mode=%s, n_init=%d, elapsed=%.3fs.",
                init_mode, kmeans_n_init, time.perf_counter() - fit_start)
    return km


def seed_kmeans_full(logger, contig_file: str, namelist: List[str], out_path: str,
                     X_mat: np.ndarray, bin_number: int, prefix: str, length_weight: np.ndarray,
                     contig_lengths: np.ndarray, seed_indexes: np.ndarray,
                     random_seed: Optional[int], minimum_length: int):
    """
    Perform weighted seed-kmeans clustering with specified parameters.

    Parameters:
    :param contig_file: The path to the contig file.
    :param namelist: A list of contig names.
    :param out_path: The output path for saving results.
    :param X_mat: The input data matrix for clustering.
    :param bin_number: The number of bins (clusters) to create.
    :param prefix: A prefix to be added to the output file names.
    :param length_weight: The weights for contig lengths.
    :param seed_indexes: Seed positions aligned to the retained contig axis.

    :return: None

    This function performs weighted seed-based k-means clustering on the input data using specified parameters and saves the results.
    """
    time_start = time.time()
    output_path = Path(out_path)
    labels_path = output_path / kmeans_label_name(bin_number, prefix)
    bins_path = output_path / kmeans_bins_name(bin_number, prefix)
    if labels_path.exists():
        labels = validate_cluster_labels(labels_path, len(namelist))
    else:
        km = fit_seed_kmeans(logger, X_mat, length_weight, bin_number, seed_indexes, random_seed)
        save_cluster_labels_atomic(labels_path, km.labels_, len(namelist))
        labels = validate_cluster_labels(labels_path, len(namelist))
    if not bins_path.exists():
        write_bins_from_labels(
            contig_file, namelist, contig_lengths, labels, bins_path,
            minimum_length=minimum_length)
    del labels
    logger.info("Running weighted seed-kmeans cost:\t"+str(time.time() - time_start) + 's.')


def load_cluster_input(emb_file):
    embedding_file = Path(emb_file).expanduser().resolve()
    if embedding_file.name != 'embeddings.npy':
        raise ValueError(f'Clustering requires the canonical embeddings.npy file: {embedding_file}')
    bundle = load_embedding_bundle(embedding_file)
    return bundle.embeddings, np.asarray(bundle.contig_ids), bundle.contig_lengths


def read_seed_names(seed_file, logger):
    seed_names = []
    seen = set()
    duplicate_count = 0

    with open(seed_file, encoding='utf-8') as input_handle:
        for line in input_handle:
            seed_name = line.rstrip('\r\n').split('\t', 1)[0]
            if not seed_name:
                continue
            if seed_name in seen:
                duplicate_count += 1
                continue
            seen.add(seed_name)
            seed_names.append(seed_name)

    if not seed_names:
        raise ValueError(f'No seed identifiers found in {seed_file}.')
    if duplicate_count:
        logger.warning("Ignored %d duplicate seed identifiers in %s.", duplicate_count, seed_file)
    return seed_names


def resolve_seed_indexes(seed_names, contig_ids, retained, logger):
    wanted = set(seed_names)
    resolved = {}
    retained_index = 0

    for source_index, contig_id in enumerate(contig_ids):
        if contig_id in wanted and contig_id not in resolved:
            resolved[contig_id] = retained_index if retained[source_index] else None
        if retained[source_index]:
            retained_index += 1

    absent = [name for name in seed_names if name not in resolved]
    filtered = [name for name in seed_names if name in resolved and resolved[name] is None]
    seed_indexes = np.asarray(
        [resolved[name] for name in seed_names if resolved.get(name) is not None], dtype=np.int64)

    if absent:
        logger.warning("%d/%d seed contigs are absent from the clustering input.", len(absent), len(seed_names))
    if filtered:
        logger.warning("%d/%d seed contigs were removed by the contig-length filter.", len(filtered), len(seed_names))
    if seed_indexes.size == 0:
        raise ValueError('No seed contigs remain on the filtered clustering axis.')
    return seed_indexes


def prepare_cluster_axis(contig_ids, contig_lengths, seed_file, contig_len, logger):
    contig_ids = np.asarray(contig_ids)
    contig_lengths = np.asarray(contig_lengths)
    if contig_ids.ndim != 1 or contig_lengths.ndim != 1 or len(contig_ids) != len(contig_lengths):
        raise ValueError('Contig identifiers and lengths must be aligned one-dimensional arrays.')

    n50 = calculate_n50(contig_lengths)
    retained = contig_lengths >= contig_len
    seed_names = read_seed_names(seed_file, logger)
    seed_indexes = resolve_seed_indexes(seed_names, contig_ids, retained, logger)
    return retained, contig_ids[retained], contig_lengths[retained], seed_indexes, len(seed_names), n50


def validate_cluster_features(features, contig_ids, contig_lengths):
    if features.ndim != 2:
        raise ValueError(f'Clustering features must be two-dimensional; received shape {features.shape}.')
    if contig_ids.ndim != 1 or contig_lengths.ndim != 1:
        raise ValueError('Contig identifiers and lengths must be one-dimensional.')
    if len(features) != len(contig_ids) or len(contig_ids) != len(contig_lengths):
        raise ValueError(
            f'Clustering arrays are not aligned: features={len(features)}, identifiers={len(contig_ids)}, '
            f'lengths={len(contig_lengths)}.')
    if len(contig_ids) == 0:
        raise ValueError('No contigs remain on the clustering axis.')


def get_kmeans_cluster_counts(seed_count, cluster_num, retained_seed_count, retained_contig_count):
    if seed_count < 1:
        raise ValueError('Seed cluster count must be greater than zero.')
    if seed_count > retained_contig_count:
        raise ValueError(
            f'Default seed cluster count {seed_count} exceeds the retained contig count '
            f'{retained_contig_count}.')

    bin_numbers = [seed_count]
    if cluster_num and cluster_num != seed_count:
        bin_numbers.append(cluster_num)

    for bin_number in bin_numbers:
        if bin_number < retained_seed_count or bin_number > retained_contig_count:
            raise ValueError(
                f'Invalid cluster count {bin_number}: retained_seeds={retained_seed_count}, '
                f'retained_contigs={retained_contig_count}.')
    return bin_numbers


def _available_cpu_count():
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def resolve_clustering_worker_counts(num_threads, leiden_workers, available_cpu_count):
    if num_threads < 1:
        raise ValueError('num_threads must be greater than zero.')
    requested_leiden_workers = num_threads if leiden_workers is None else leiden_workers
    if requested_leiden_workers < 1:
        raise ValueError('leiden_workers must be greater than zero.')
    return min(num_threads, available_cpu_count), requested_leiden_workers, \
        min(requested_leiden_workers, available_cpu_count)


def validate_max_edges(max_edges, contig_count):
    if max_edges < 1 or max_edges >= contig_count:
        raise ValueError(
            f'max_edges must be between 1 and the retained contig count minus one; '
            f'received max_edges={max_edges}, retained_contigs={contig_count}.')


def _read_proc_value_mib(path, key):
    try:
        with open(path) as input_handle:
            for line in input_handle:
                name, _, value = line.partition(':')
                if name == key:
                    return int(value.split()[0]) / 1024
    except (OSError, ValueError, IndexError):
        pass
    return float('nan')


def _process_memory_mib():
    return _read_proc_value_mib('/proc/self/status', 'VmRSS'), _read_proc_value_mib('/proc/self/smaps_rollup', 'Pss')


def _available_memory_mib():
    return _read_proc_value_mib('/proc/meminfo', 'MemAvailable')


def prepare_leiden_context(namelist, ann_neighbor_indices, ann_distances, length_weight,
                           max_edges, is_membership_fixed, jobs, lmode='l2'):
    """Prepare all read-only Leiden graph families before forking one worker pool."""
    vertex_count = len(namelist)
    sources = np.repeat(np.arange(vertex_count, dtype=np.int64), max_edges)
    targets = ann_neighbor_indices[:, 1:].reshape(-1)
    distances = ann_distances[:, 1:].reshape(-1)

    if len(sources) != len(targets) or len(targets) != len(distances):
        raise ValueError('HNSW neighbor arrays do not match max_edges.')

    if lmode not in ('l1', 'l2'):
        raise ValueError('Unsupported Leiden distance mode: ' + str(lmode))

    required_bandwidths = {}
    for _, partgraph_ratio, bandwidth, _ in jobs:
        ratio_bandwidths = required_bandwidths.setdefault(partgraph_ratio, [])
        if bandwidth not in ratio_bandwidths:
            ratio_bandwidths.append(bandwidth)

    graphs = {}
    weights = {}
    for partgraph_ratio, bandwidths in required_bandwidths.items():
        dist_cutoff = np.percentile(distances, partgraph_ratio)
        save_index = (distances <= dist_cutoff) & (sources > targets)
        edge_distances = distances[save_index]
        edges = np.column_stack((sources[save_index], targets[save_index]))
        graphs[partgraph_ratio] = Graph(vertex_count, edges, directed=False)

        if lmode == 'l1':
            edge_distances = np.sqrt(edge_distances)

        for bandwidth in bandwidths:
            edge_weights = np.exp(-edge_distances / bandwidth)
            edge_weights.setflags(write=False)
            weights[(partgraph_ratio, bandwidth)] = edge_weights

    return {
        'expected_rows': vertex_count,
        'length_weight': length_weight,
        'is_membership_fixed': is_membership_fixed,
        'graphs': graphs,
        'weights': weights,
    }


def run_leiden(job):
    """Run one Leiden parameter combination using inputs inherited through fork."""
    if _LEIDEN_CONTEXT is None:
        raise RuntimeError('Leiden worker context is not initialized.')

    output_file, partgraph_ratio, bandwidth, resolution_parameter = job
    graph = _LEIDEN_CONTEXT['graphs'][partgraph_ratio]
    edge_weights = _LEIDEN_CONTEXT['weights'][(partgraph_ratio, bandwidth)]
    res = leidenalg.RBERVertexPartition(
        graph, weights=edge_weights, initial_membership=None,
        resolution_parameter=resolution_parameter, node_sizes=_LEIDEN_CONTEXT['length_weight'])

    optimiser = leidenalg.Optimiser()
    optimiser.consider_comms = leidenalg.ALL_NEIGH_COMMS
    optimiser.refine_consider_comms = leidenalg.ALL_NEIGH_COMMS
    if _LEIDEN_CONTEXT['seed'] is not None:
        optimiser.set_rng_seed(_LEIDEN_CONTEXT['seed'])
    time_start = time.monotonic()
    total_improvement = optimiser.optimise_partition(
        res, is_membership_fixed=_LEIDEN_CONTEXT['is_membership_fixed'], n_iterations=-1)

    rows = save_cluster_labels_atomic(
        output_file, res.membership, _LEIDEN_CONTEXT['expected_rows'])

    worker_rss_mib, worker_pss_mib = _process_memory_mib()
    return output_file, partgraph_ratio, bandwidth, resolution_parameter, total_improvement, \
        time.monotonic() - time_start, rows, os.getpid(), worker_rss_mib, worker_pss_mib


def run_leiden_jobs(logger, context, jobs, num_workers):
    """Run all pending Leiden job tuples in one persistent fork pool."""
    global _LEIDEN_CONTEXT

    if not jobs:
        return []
    if num_workers < 1:
        raise ValueError('Leiden worker count must be greater than zero.')

    effective_workers = min(num_workers, len(jobs))
    fork_context = multiprocessing.get_context('fork')
    _LEIDEN_CONTEXT = context
    multiprocess = None
    results = []

    try:
        multiprocess = fork_context.Pool(effective_workers)
        for result in multiprocess.imap_unordered(run_leiden, jobs, chunksize=1):
            results.append(result)
            output_file, ratio, bandwidth, resolution, improvement, elapsed, rows, worker_pid, \
                worker_rss_mib, worker_pss_mib = result
            logger.info(
                "Leiden completed %d/%d: ratio=%s, bandwidth=%s, resolution=%s, improvement=%.12g, elapsed=%.3fs, rows=%d, worker_pid=%d, worker_rss_mib=%.1f, worker_pss_mib=%.1f, output=%s.",
                len(results), len(jobs), ratio, bandwidth, resolution, improvement, elapsed, rows,
                worker_pid, worker_rss_mib, worker_pss_mib, output_file)
        multiprocess.close()
        multiprocess.join()
    except BaseException:
        if multiprocess is not None:
            multiprocess.terminate()
            multiprocess.join()
        raise
    finally:
        _LEIDEN_CONTEXT = None

    return results


def cluster(logger, args, prefix=None):
    """Load the canonical clustering input and dispatch the aligned array core."""
    logger.info("Start clustering.")
    features, contig_ids, contig_lengths = load_cluster_input(args.emb_file)
    logger.info("Loaded clustering input: mode=binary-bundle, contigs=%d, features=%d, dtype=%s.",
                len(contig_ids), features.shape[1], features.dtype)
    retained, contig_ids, contig_lengths, seed_indexes, seed_count, n50 = prepare_cluster_axis(
        contig_ids, contig_lengths, args.seed_file, args.contig_len, logger)
    logger.info(
        "Prepared clustering axis: retained_contigs=%d, unique_seeds=%d, retained_seeds=%d.",
        len(contig_ids), seed_count, len(seed_indexes))
    if not retained.all():
        features = features[retained]
    return cluster_features(
        logger, args, features, contig_ids, contig_lengths, seed_indexes, seed_count, n50, prefix)


def cluster_features(logger, args, features, contig_ids, contig_lengths,
                     seed_indexes, seed_count, n50, prefix=None):
    """Cluster complete aligned arrays without reading feature, FASTA, or seed files."""
    features = np.asarray(features)
    contig_ids = np.asarray(contig_ids)
    contig_lengths = np.asarray(contig_lengths)
    seed_indexes = np.asarray(seed_indexes, dtype=np.int64)
    validate_cluster_features(features, contig_ids, contig_lengths)
    if seed_indexes.ndim != 1 or seed_indexes.size == 0:
        raise ValueError('Retained seed indexes must be a nonempty one-dimensional array.')
    if seed_indexes.min() < 0 or seed_indexes.max() >= len(contig_ids):
        raise ValueError('Retained seed indexes fall outside the clustering axis.')

    logger.info('N50:\t' + str(n50))
    length_weight = contig_lengths.tolist()
    output_path = args.output_path + '/cluster_res/'

    if args.not_l2normaize:
        norm_embeddings = features
    else:
        norm_embeddings = normalize(features)
    os.makedirs(output_path, exist_ok=True)

    mode = 'weight_seed_kmeans'
    if prefix:
        mode = mode + '_' + prefix
    logger.info("Run weighted seed k-means for obtaining the SCG information of the contigs within a manageable time during the final step.")
    bin_nums = get_kmeans_cluster_counts(seed_count, args.cluster_num, len(seed_indexes), len(contig_ids))

    logger.info("Bin_numbers:\t"+str(bin_nums))
    for bin_number in bin_nums:
        logger.info(bin_number)
        seed_kmeans_full(logger, args.contig_file, contig_ids, output_path, norm_embeddings,
                         bin_number, mode, length_weight, contig_lengths, seed_indexes,
                         args.seed, args.contig_len)

    available_cpu_count = _available_cpu_count()
    hnsw_num_threads, requested_leiden_workers, leiden_num_workers = resolve_clustering_worker_counts(
        args.num_threads, args.leiden_workers, available_cpu_count)
    max_edges = args.max_edges
    validate_max_edges(max_edges, len(contig_ids))
    logger.info(
        "Clustering workers: requested_hnsw_threads=%d, effective_hnsw_threads=%d, requested_leiden_workers=%d, effective_leiden_workers=%d, available_cpus=%d.",
        args.num_threads, hnsw_num_threads, requested_leiden_workers, leiden_num_workers, available_cpu_count)
    if args.seed is not None:
        logger.info("Clustering reproducibility: seed=%d.", args.seed)
    else:
        logger.info("Clustering reproducibility: random seed is not fixed.")
    logger.info(
        "Leiden settings: consider_comms=ALL_NEIGH_COMMS, refine_consider_comms=ALL_NEIGH_COMMS, n_iterations=-1 (until no improvement).")

    jobs = []
    for partgraph_ratio in LEIDEN_PARTGRAPH_RATIOS:
        for bandwidth in LEIDEN_BANDWIDTHS:
            for resolution in LEIDEN_RESOLUTIONS:
                output_file = os.path.join(
                    output_path, leiden_output_name(max_edges, partgraph_ratio, bandwidth, resolution))
                if os.path.exists(output_file):
                    labels = validate_cluster_labels(output_file, len(contig_ids))
                    del labels
                else:
                    jobs.append((output_file, partgraph_ratio, bandwidth, resolution))

    if jobs:
        p = fit_hnsw_index(logger, norm_embeddings, num_threads=hnsw_num_threads,
                           ef=max_edges * 10, random_seed=args.seed)
        time_start = time.time()
        logger.info("Query completed HNSW index with %d threads.", hnsw_num_threads)
        ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges + 1, num_threads=hnsw_num_threads)
        time_end = time.time()
        logger.info('knn query time cost:\t' + str(time_end - time_start) + "s")
        del p

        is_membership_fixed = np.zeros(len(contig_ids), dtype=bool)
        is_membership_fixed[seed_indexes] = True
        is_membership_fixed = is_membership_fixed.tolist()

        context = prepare_leiden_context(
            contig_ids, ann_neighbor_indices, ann_distances, length_weight, max_edges,
            is_membership_fixed, jobs, lmode='l2')
        del ann_neighbor_indices, ann_distances
        context['seed'] = args.seed
        for partgraph_ratio, graph in context['graphs'].items():
            logger.info(
                "Leiden graph resources: ratio=%s, vertices=%d, edges=%d.",
                partgraph_ratio, graph.vcount(), graph.ecount())
        effective_pool_workers = min(leiden_num_workers, len(jobs))
        parent_rss_mib, parent_pss_mib = _process_memory_mib()
        logger.info(
            "Leiden pool resources: jobs=%d, requested_workers=%d, effective_workers=%d, graph_families=%d, weight_families=%d, parent_rss_mib=%.1f, parent_pss_mib=%.1f, available_memory_mib=%.1f.",
            len(jobs), requested_leiden_workers, effective_pool_workers, len(context['graphs']),
            len(context['weights']), parent_rss_mib, parent_pss_mib, _available_memory_mib())
        pool_start = time.monotonic()
        results = run_leiden_jobs(logger, context, jobs, effective_pool_workers)
        logger.info(
            "Leiden pool completed: completed=%d/%d, elapsed=%.3fs.",
            len(results), len(jobs), time.monotonic() - pool_start)
        del results, context
        gc.collect()
        released_rss_mib, released_pss_mib = _process_memory_mib()
        logger.info(
            "Leiden pool released: parent_rss_mib=%.1f, parent_pss_mib=%.1f, available_memory_mib=%.1f.",
            released_rss_mib, released_pss_mib, _available_memory_mib())
    else:
        logger.info('All Leiden outputs already exist for max_edges=' + str(max_edges) + '.')

    validated_names, validated_rows = validate_leiden_output_set(
        output_path, max_edges, expected_rows=len(contig_ids))
    logger.info(
        "Validated complete Leiden result grid: max_edges=%d, files=%d, rows_per_file=%d.",
        max_edges, len(validated_names), validated_rows)
