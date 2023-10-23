import numpy as np
import pandas as pd
import torch
from utils import get_kmerMetric_emb
from sklearn.preprocessing import normalize


def get_kmer_coverage(data_path: str, n_views: int = 2, kmer_model_path: str = 'empty',
                      device = torch.device('cpu'), nokmer: bool = False, cov_meannormalize: bool = False,
                      cov_minmaxnormalize: bool = False, cov_standardization: bool = False, addvars: bool = False,
                      vars_sqrt: bool = False, kmer_l2_normalize: bool = False, kmerMetric_notl2normalize: bool = False):
    """
    Get features

    :param data_path: The path to the data directory.
    :param n_views: The number of views (default: 2).
    :param kmer_model_path: The path to the k-mer model (default: 'empty').
    :param device: The device for computation (default: 'cpu').
    :param nokmer: Flag to exclude k-mer data (default: False).
    :param cov_meannormalize: Flag to mean normalize coverage data (default: False).
    :param cov_minmaxnormalize: Flag to min-max normalize coverage data (default: False).
    :param cov_standardization: Flag to standardize coverage data (default: False).
    :param addvars: Flag to include additional variables (default: False).
    :param vars_sqrt: Flag to take the square root of variables (default: False).
    :param kmer_l2_normalize: Flag to L2 normalize k-mer data (default: False).
    :param kmerMetric_notl2normalize: Flag to not L2 normalize k-mer metric data (default: False).

    :return: A list of preprocessed data and a list of contig names.
    """
    namelist = pd.read_csv(data_path + 'aug0_datacoverage_mean.tsv', sep='\t', usecols=range(1)).values[:, 0]

    mapObj = dict(zip(namelist, range(len(namelist))))
    for view in range(n_views):
        cov_file = data_path + 'aug' + str(view) + '_datacoverage_mean.tsv'
        if not nokmer:
            com_file = data_path + 'aug' + str(view) + '/kmer_4_f0.csv'

        covHeader = pd.read_csv(cov_file, sep='\t', nrows=1)
        shuffled_covMat = pd.read_csv(cov_file, sep='\t', usecols=range(1, covHeader.shape[1])).values
        shuffled_namelist = pd.read_csv(cov_file, sep='\t', usecols=range(1)).values[:, 0]

        covIdxArr = np.empty(len(mapObj), dtype=np.int)
        for contigIdx in range(len(shuffled_namelist)):
            if shuffled_namelist[contigIdx].split('_aug')[0] in mapObj:
                covIdxArr[mapObj[shuffled_namelist[contigIdx].split('_aug')[0]]] = contigIdx
        covMat = shuffled_covMat[covIdxArr]


        if not nokmer:
            compositHeader = pd.read_csv(com_file, sep=',', nrows=1)
            shuffled_compositMat = pd.read_csv(com_file, sep=',', usecols=range(1, compositHeader.shape[1])).values
            shuffled_namelist = pd.read_csv(com_file, sep=',', usecols=range(1)).values[:, 0]

            covIdxArr = np.empty(len(mapObj), dtype=np.int)
            for contigIdx in range(len(shuffled_namelist)):
                if shuffled_namelist[contigIdx].split('_aug')[0] in mapObj:
                    covIdxArr[mapObj[shuffled_namelist[contigIdx].split('_aug')[0]]] = contigIdx
            compositMat = shuffled_compositMat[covIdxArr]

        if addvars:
            vars_file = data_path+'aug'+str(view)+'_datacoverage_var.tsv'

            varsHeader = pd.read_csv(vars_file, sep='\t', nrows=1)
            shuffled_varsMat = pd.read_csv(vars_file, sep='\t', usecols=range(1, varsHeader.shape[1])).values
            shuffled_namelist = pd.read_csv(vars_file, sep='\t', usecols=range(1)).values[:, 0]

            covIdxArr = np.empty(len(mapObj), dtype=np.int)
            for contigIdx in range(len(shuffled_namelist)):
                if shuffled_namelist[contigIdx].split('_aug')[0] in mapObj:
                    covIdxArr[mapObj[shuffled_namelist[contigIdx].split('_aug')[0]]] = contigIdx
            varsMat = shuffled_varsMat[covIdxArr]

        if view == 0:
            covMats = covMat
            if addvars:
                varsMats = varsMat
            if not nokmer:
                compositMats = compositMat

        else:
            covMats = np.vstack((covMats, covMat))
            if addvars:
                varsMats = np.vstack((varsMats, varsMat))
            if not nokmer:
                compositMats = np.vstack((compositMats, compositMat))

    # use cov_maxnormalize
    covMats = covMats + 1e-5
    if addvars:
        if vars_sqrt:
            varsMats = np.sqrt(varsMats)
        varsMats = varsMats + 1e-5

    if cov_meannormalize:
        covMats = covMats / covMats.mean(axis=0)[None, :]
        if addvars:
            varsMats = varsMats / varsMats.mean(axis=0)[None, :]

    elif cov_minmaxnormalize:
        covMats = (covMats - covMats.min(axis=0)[None, :])/ (covMats.max(axis=0)[None, :] - covMats.min(axis=0)[None, :] )
        if addvars:
            varsMats = (varsMats - varsMats.min(axis=0)[None, :])/ (varsMats.max(axis=0)[None, :] - varsMats.min(axis=0)[None, :] )

    elif cov_standardization:
        covMats = (covMats - covMats.mean(axis=0)[None, :])/ covMats.std(axis=0)[None, :]
        if addvars:
            varsMats = (varsMats - varsMats.mean(axis=0)[None, :])/ varsMats.std(axis=0)[None, :]
    else:
        covMats = covMats / covMats.max(axis=0)[None, :]
        if addvars:
            varsMats = varsMats / varsMats.max(axis=0)[None, :]

    if not nokmer:
        compositMats = compositMats + 1
        compositMats = compositMats / compositMats.sum(axis=1)[:, None]
        if kmer_l2_normalize:
            compositMats = normalize(compositMats)
        if kmer_model_path != 'empty':
            compositMats = get_kmerMetric_emb(kmer_model_path, compositMats, device,kmerMetric_notl2normalize)

    if not nokmer:
        X_ts = np.hstack((covMats, compositMats))
    else:
        X_ts = covMats

    if addvars:
        X_ts = np.hstack((varsMats, X_ts))

    return list(torch.split(torch.from_numpy(X_ts).float(), len(namelist))), namelist


def get_ContrastiveLearningDataset(data_path: str, n_views: int = 2, kmer_model_path: str = 'empty',
                                   device=torch.device('cpu'), nokmer: bool = False, cov_meannormalize: bool = False,
                                   cov_minmaxnormalize: bool = False, cov_standardization: bool = False, addvars: bool = False,
                                   vars_sqrt: bool = False, kmer_l2_normalize: bool = False, kmerMetric_notl2normalize: bool = False):
    """
    Get a Contrastive Learning dataset based on input parameters.

    :param data_path: The path to the data.
    :param n_views: The number of views for data (default: 2).
    :param kmer_model_path: The path to the k-mer model (default: 'empty').
    :param device: The device to use for computations (default: 'cpu').
    :param nokmer: Whether to use k-mer features (default: False).
    :param cov_meannormalize: Whether to mean normalize coverage (default: False).
    :param cov_minmaxnormalize: Whether to min-max normalize coverage (default: False).
    :param cov_standardization: Whether to standardize coverage (default: False).
    :param addvars: Whether to add additional variables (default: False).
    :param vars_sqrt: Whether to take the square root of additional variables (default: False).
    :param kmer_l2_normalize: Whether to L2 normalize k-mer features (default: False).
    :param kmerMetric_notl2normalize: Whether not to L2 normalize k-mer features (default: False).

    :return: A tuple containing the dataset and a list of names.
    """
    if not data_path.endswith('/'):
        data_path = data_path + '/'
    dataset, namelist = get_kmer_coverage(data_path, n_views, kmer_model_path, device, nokmer, cov_meannormalize, cov_minmaxnormalize, cov_standardization,addvars,vars_sqrt,kmer_l2_normalize,kmerMetric_notl2normalize)
    return dataset, namelist
