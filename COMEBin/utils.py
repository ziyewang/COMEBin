import gzip
from Bio import SeqIO
import mimetypes
import os
import shutil
import pandas as pd
import numpy as np

import torch
import yaml
import sys

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


def get_kmer_coverage_aug0(data_path):
    namelist = pd.read_csv(data_path + 'aug0_datacoverage_mean.tsv', sep='\t', usecols=range(1)).values[:, 0]
    mapObj = dict(zip(namelist, range(len(namelist))))

    cov_file = data_path + 'aug0_datacoverage_mean.tsv'
    com_file = data_path + 'aug0/kmer_4_f0.csv'

    covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
    shuffled_covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).values
    shuffled_namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]

    covIdxArr = np.empty(len(mapObj), dtype=np.int)
    for contigIdx in range(len(shuffled_namelist)):
        if shuffled_namelist[contigIdx].split('_aug')[0] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx].split('_aug')[0]]] = contigIdx
    covMat = shuffled_covMat[covIdxArr]

    compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
    shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).values
    shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]

    covIdxArr = np.empty(len(mapObj), dtype=np.int)
    for contigIdx in range(len(shuffled_namelist)):
        if shuffled_namelist[contigIdx].split('_aug')[0] in mapObj:
            covIdxArr[mapObj[shuffled_namelist[contigIdx].split('_aug')[0]]] = contigIdx
    compositMat = shuffled_compositMat[covIdxArr]

    # use cov_maxnormalize
    covMat = covMat + 1e-5
    covMat = covMat / covMat.max(axis=0)[None, :]

    compositMat = compositMat + 1
    compositMat = compositMat / compositMat.sum(axis=1)[:, None]

    X_t = np.hstack((covMat, compositMat))

    return X_t, covMat, compositMat, namelist




@torch.no_grad()
def get_embeddings(model, test_dl, embedding_dim, device):
    model.eval()
    embs = np.zeros(shape=(0, embedding_dim))
    for x in test_dl:
        emb = model(x.to(device)).to('cpu')
        embs = np.vstack((embs, emb))

    return embs


def get_kmerMetric_emb(kmer_model_path,compositMats,device=torch.device('cpu'),kmerMetric_notl2normalize=False):
    #####
    # load kmer_metric
    config_file = os.path.dirname(kmer_model_path) + '/kmerMetric_config.yaml'

    from ruamel.yaml import YAML
    from pathlib import Path
    import torch.nn as nn
    from models.mlp import EmbeddingNet
    from sklearn.preprocessing import normalize


    yaml = YAML(typ='safe')

    cnf = yaml.load(Path(config_file))

    ps = [cnf['dropout_value']] * (len(cnf['emb_szs']) - 1)
    actn = nn.LeakyReLU()

    kmerMetric_model = EmbeddingNet(
        in_sz=len(compositMats[0]),
        out_sz=cnf['embedding_dim'],
        emb_szs=cnf['emb_szs'],
        ps=ps,
        use_bn=True,
        actn=actn,
    )

    kmerMetric_model = kmerMetric_model.to(device)

    kmerMetric_model.load_state_dict(torch.load(kmer_model_path, map_location=device))

    compositMats = torch.from_numpy(compositMats).float()

    test_dataloader = torch.utils.data.DataLoader(compositMats,
                                                  batch_size=2048,
                                                  shuffle=False)

    compositMats = get_embeddings(kmerMetric_model, test_dataloader, cnf['embedding_dim'], device)
    if not kmerMetric_notl2normalize:
        compositMats = normalize(compositMats)
    return compositMats


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def gen_seed(logger, contig_file: str, threads: int, contig_length_threshold: int,
             marker_name: str = "marker", quarter: str = "3quarter"):
    """
    Generate seed sequences from contigs using FragGeneScan, HMMsearch, and custom markers.

    :param contig_file: Path to the input contig file.
    :param threads: The number of threads to use for processing.
    :param contig_length_threshold: The contig length threshold.
    :param marker_name: The marker name (default: "marker").
    :param quarter: The quarter identifier (default: "3quarter").
    :return: The number of candidate seeds generated.
    """
    fragScanURL = 'run_FragGeneScan.pl'

    hmmExeURL = 'hmmsearch'
    markerExeURL = os.path.join(os.getcwd(), '../auxiliary', 'test_getmarker_' + quarter + '.pl')
    markerURL = os.path.join(os.getcwd(), '../auxiliary', marker_name + '.hmm')
    seedURL = contig_file + "." + marker_name + "." + quarter + "_lencutoff_" + str(contig_length_threshold) + ".seed"
    fragResultURL = contig_file + ".frag.faa"
    hmmResultURL = contig_file + '.' + marker_name + ".hmmout"

    if not (os.path.exists(fragResultURL)):
        fragCmd = fragScanURL + " -genome=" + contig_file + " -out=" + contig_file + ".frag -complete=0 -train=complete -thread=" + str(
            threads) + " 1>" + contig_file + ".frag.out 2>" + contig_file + ".frag.err"
        logger.info("exec cmd: " + fragCmd)
        os.system(fragCmd)

    if os.path.exists(fragResultURL):
        if not (os.path.exists(hmmResultURL)):
            hmmCmd = hmmExeURL + " --domtblout " + hmmResultURL + " --cut_tc --cpu " + str(
                threads) + " " + markerURL + " " + fragResultURL + " 1>" + hmmResultURL + ".out 2>" + hmmResultURL + ".err"
            logger.info("exec cmd: " + hmmCmd)
            os.system(hmmCmd)

        if os.path.exists(hmmResultURL):
            if not (os.path.exists(seedURL)):
                markerCmd = markerExeURL + " " + hmmResultURL + " " + contig_file + " " + str(
                    contig_length_threshold) + " " + seedURL
                logger.info("exec cmd: " + markerCmd)
                os.system(markerCmd)

            if os.path.exists(seedURL):
                candK = file_len(seedURL)
            else:
                logger.info("markerCmd failed! Not exist: " + markerCmd)
                candK = 0
        else:
            logger.info("Hmmsearch failed! Not exist: " + hmmResultURL)
            sys.exit()
    else:
        logger.info("FragGeneScan failed! Not exist: " + fragResultURL)
        sys.exit()
    return candK