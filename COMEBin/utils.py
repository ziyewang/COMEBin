import gzip
from Bio import SeqIO
from collections import deque
import mimetypes
import os
from pathlib import Path
import random
import shlex
import shutil
import subprocess
import tempfile
from typing import Optional
import numpy as np

import torch
import yaml


DEFAULT_HMM_EVALUE = 1e-5
DEFAULT_RANDOM_SEED = None
COMEBIN_DIR = Path(__file__).resolve().parent
AUXILIARY_DIR = COMEBIN_DIR.parent / 'auxiliary'


def _process_memory_status():
    status = {}
    try:
        with open('/proc/self/status') as status_file:
            for line in status_file:
                key, _, value = line.partition(':')
                if key in ('VmRSS', 'VmPeak'):
                    status[key] = value.strip()
    except OSError:
        pass
    return status


def configure_reproducibility(seed: Optional[int]) -> None:
    if seed is None:
        return
    if seed < 0 or seed >= 2**32:
        raise ValueError('seed must be between 0 and 2**32 - 1.')

    os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)


def seed_data_loader_worker(worker_id: int) -> None:
    worker_seed = torch.initial_seed() % 2**32
    random.seed(worker_seed)
    np.random.seed(worker_seed)


def save_checkpoint(state, is_best, filename='checkpoint.pth.tar'):
    torch.save(state, filename)
    if is_best:
        shutil.copyfile(filename, 'model_best.pth.tar')


def save_config_file(model_checkpoints_folder, args):
    if not os.path.exists(model_checkpoints_folder):
        os.makedirs(model_checkpoints_folder)
    with open(os.path.join(model_checkpoints_folder, 'config.yml'), 'w') as outfile:
        yaml.dump(args, outfile, default_flow_style=False)


def accuracy(output, target, topk=(1,)):
    """Computes the accuracy over the k top predictions for the specified values of k"""
    with torch.no_grad():
        maxk = max(topk)
        batch_size = target.size(0)

        _, pred = output.topk(maxk, 1, True, True)
        pred = pred.t()
        correct = pred.eq(target.view(1, -1).expand_as(pred))

        res = []
        for k in topk:
            correct_k = correct[:k].reshape(-1).float().sum(0, keepdim=True)
            res.append(correct_k.mul_(100.0 / batch_size))
        return res


def get_length(fastx_file):
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))
    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length
    file_format = None
    line = f.readline()
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
    f.seek(0)
    if not file_format:
        raise RuntimeError("Invalid sequence file: '{}".format(fastx_file))
    for seq_record in SeqIO.parse(f, file_format):
        length[seq_record.id] = len(seq_record.seq)

    f.close()
    return length



def save_result(result, filepath, namelist):
    filedir, filename = os.path.split(filepath)
    if not filename:
        filename = "result.tsv"
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    f = open(filepath, 'w')
    for contigIdx in range(len(result)):
        f.write(namelist[contigIdx] + "\t" + str(result[contigIdx].item(0)) + "\n")
    f.close()


def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    N50 = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break
    return N50


def _require_file(path, label):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")
    return path


def _require_executable(name):
    executable = shutil.which(name)
    if executable is None:
        raise FileNotFoundError(f"Required executable is not available on PATH: {name}")
    return executable


def _require_sibling_executable(executable, sibling_name):
    sibling = Path(executable).parent / sibling_name
    if not sibling.is_file() or not os.access(sibling, os.X_OK):
        raise FileNotFoundError(f"Required executable is not available beside {executable}: {sibling_name}")
    return str(sibling.resolve())


def _stderr_tail(path, line_count=20):
    try:
        with open(path) as input_handle:
            return ''.join(deque(input_handle, maxlen=line_count)).strip()
    except OSError:
        return ''


def _run_external_command(logger, label, command, stdout_path, stderr_path):
    command = [str(argument) for argument in command]
    logger.info("%s command: %s", label, shlex.join(command))
    try:
        with open(stdout_path, 'w') as stdout_handle, open(stderr_path, 'w') as stderr_handle:
            completed = subprocess.run(command, stdout=stdout_handle, stderr=stderr_handle, check=False)
    except OSError as error:
        raise RuntimeError(f"{label} could not be started: {error}") from error
    logger.info("%s exited with status %d; stderr=%s.", label, completed.returncode, stderr_path)
    if completed.returncode != 0:
        stderr_tail = _stderr_tail(stderr_path)
        detail = f"; tail={stderr_tail}" if stderr_tail else ''
        raise RuntimeError(f"{label} failed with status {completed.returncode}; stderr={stderr_path}{detail}")
    return completed


def _temporary_result_path(final_path):
    final_path = Path(final_path)
    descriptor, temp_name = tempfile.mkstemp(prefix=f'.{final_path.name}.{os.getpid()}.', suffix='.tmp',
                                             dir=final_path.parent)
    os.close(descriptor)
    return Path(temp_name)


def _count_hmm_hits(hmm_path):
    with open(hmm_path) as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith('#'))


def _count_seed_lines(seed_path):
    with open(seed_path) as handle:
        return sum(1 for line in handle if line.strip())


def _hmm_models_have_tc_cutoffs(hmm_file):
    model_count = 0
    tc_count = 0
    with open(hmm_file) as handle:
        for line in handle:
            if line.startswith('NAME'):
                model_count += 1
            elif line.startswith('TC'):
                tc_count += 1
    return model_count > 0 and tc_count == model_count


def gen_seed(logger, contig_file: str, threads: int, contig_length_threshold: int,
             marker_name: str = "marker", quarter: str = "3quarter",
             hmm_evalue: float = DEFAULT_HMM_EVALUE):
    """
    Generate seed sequences from contigs using FragGeneScan, HMMsearch, and custom markers.

    :param contig_file: Path to the input contig file.
    :param threads: The number of threads to use for processing.
    :param contig_length_threshold: The contig length threshold.
    :param marker_name: The marker name (default: "marker").
    :param quarter: The quarter identifier (default: "3quarter").
    :param hmm_evalue: Sequence E-value used when the marker HMM does not define TC cutoffs.
    :return: The number of candidate seeds generated.
    """
    hmm_evalue = float(hmm_evalue)
    if not np.isfinite(hmm_evalue) or hmm_evalue <= 0:
        raise ValueError(f"hmm_evalue must be a finite positive number, got {hmm_evalue}.")

    seed_result = Path(contig_file + f'.{marker_name}.{quarter}_lencutoff_{contig_length_threshold}.seed')
    if seed_result.is_file() and seed_result.stat().st_size > 0:
        seed_count = _count_seed_lines(seed_result)
        if seed_count == 0:
            raise RuntimeError(f"Existing marker seed file contains no usable seeds: {seed_result}")
        logger.info("Reusing marker seed file %s with %d seeds; --hmm_evalue applies only when HMMsearch runs and cached thresholds are not compared.",
                    seed_result, seed_count)
        return seed_count

    fragscan_wrapper = _require_executable('run_FragGeneScan.pl')
    # The upstream wrapper rebuilds a shell command, so invoke its sibling binary directly to preserve path arguments.
    fragscan_executable = _require_sibling_executable(fragscan_wrapper, 'FragGeneScan')
    hmmsearch_executable = _require_executable('hmmsearch')
    perl_executable = _require_executable('perl')
    marker_script = _require_file(AUXILIARY_DIR / f'test_getmarker_{quarter}.pl', 'Marker script')
    marker_hmm = _require_file(AUXILIARY_DIR / f'{marker_name}.hmm', 'Marker HMM database')
    logger.info("Marker tools: FragGeneScan=%s, wrapper=%s, hmmsearch=%s, perl=%s, marker_script=%s, marker_hmm=%s.",
                fragscan_executable, fragscan_wrapper, hmmsearch_executable, perl_executable, marker_script, marker_hmm)

    frag_prefix = contig_file + '.frag'
    frag_result = Path(frag_prefix + '.faa')
    if not frag_result.is_file() or frag_result.stat().st_size == 0:
        _run_external_command(
            logger, 'FragGeneScan',
            [fragscan_executable, '-s', contig_file, '-o', frag_prefix, '-w', '0', '-t', 'complete', '-p', threads],
            contig_file + '.frag.out', contig_file + '.frag.err')
        if not frag_result.is_file() or frag_result.stat().st_size == 0:
            raise RuntimeError(f"FragGeneScan completed but produced no proteins: {frag_result}")
    else:
        logger.info("Reusing FragGeneScan protein file: %s.", frag_result)

    hmm_result = Path(contig_file + f'.{marker_name}.hmmout')
    if hmm_result.is_file() and hmm_result.stat().st_size > 0:
        hit_count = _count_hmm_hits(hmm_result)
        if hit_count == 0:
            raise RuntimeError(f"Existing HMMsearch result contains zero marker hits: {hmm_result}")
        logger.info("Reusing HMMsearch result %s with %d hits; --hmm_evalue applies only when HMMsearch runs and cached thresholds are not compared.",
                    hmm_result, hit_count)
    else:
        if _hmm_models_have_tc_cutoffs(marker_hmm):
            cutoff_arguments = ['--cut_tc']
            logger.info("Marker HMM %s defines TC cutoffs for every model; using --cut_tc and ignoring --hmm_evalue=%g.",
                        marker_hmm, hmm_evalue)
        else:
            hmm_evalue_argument = format(hmm_evalue, '.15g')
            cutoff_arguments = ['-E', hmm_evalue_argument]
            logger.warning("Marker HMM %s does not define TC cutoffs for every model; using -E %s.",
                           marker_hmm, hmm_evalue_argument)

        hmm_temp = _temporary_result_path(hmm_result)
        try:
            _run_external_command(
                logger, 'HMMsearch',
                [hmmsearch_executable, '--domtblout', hmm_temp, *cutoff_arguments, '--cpu', threads,
                 marker_hmm, frag_result],
                str(hmm_result) + '.out', str(hmm_result) + '.err')
            hit_count = _count_hmm_hits(hmm_temp)
            if hit_count == 0:
                raise RuntimeError(f"HMMsearch completed but found zero marker hits: {hmm_temp}")
            os.replace(hmm_temp, hmm_result)
            logger.info("Published HMMsearch result %s with %d hits.", hmm_result, hit_count)
        finally:
            if hmm_temp.exists():
                hmm_temp.unlink()

    seed_temp = _temporary_result_path(seed_result)
    try:
        _run_external_command(
            logger, 'Marker seed generation',
            [perl_executable, marker_script, hmm_result, contig_file, contig_length_threshold, seed_temp],
            str(seed_result) + '.out', str(seed_result) + '.err')
        seed_count = _count_seed_lines(seed_temp)
        if seed_count == 0:
            raise RuntimeError(f"Marker parsing completed but produced an empty seed set; hmm_hits={hit_count}; stderr={seed_result}.err")
        os.replace(seed_temp, seed_result)
        logger.info("Published marker seed file %s with %d seeds.", seed_result, seed_count)
    finally:
        if seed_temp.exists():
            seed_temp.unlink()
    return seed_count
