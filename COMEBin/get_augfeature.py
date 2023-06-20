import numpy as np
import pandas as pd
import torch
from utils import get_kmerMetric_emb
from sklearn.preprocessing import normalize


def get_kmer_coverage(data_path, kmer='4mer', n_views=2, kmer_model_path='empty',
                      device=torch.device('cpu'), nokmer=False, cov_meannormalize=False, cov_minmaxnormalize=False, cov_standardization=False,addvars=False,vars_sqrt=False,kmer_l2_normalize=False,kmerMetric_notl2normalize=False):
    # namelist = pd.read_csv(data_path + 'aug0_datacoverage_mean_edge75.tsv', sep='\t', usecols=range(1)).values[:, 0]
    namelist = pd.read_csv(data_path + 'aug0_datacoverage_mean.tsv', sep='\t', usecols=range(1)).values[:, 0]

    mapObj = dict(zip(namelist, range(len(namelist))))
    for view in range(n_views):
        # cov_file = data_path + 'aug' + str(view) + '_datacoverage_mean_edge75.tsv'
        cov_file = data_path + 'aug' + str(view) + '_datacoverage_mean.tsv'
        if not nokmer:
            if kmer == '4mer':
                com_file = data_path + 'aug' + str(view) + '/kmer_4_f0.csv'
            elif kmer == '345mer':
                com_file = data_path + 'aug' + str(view) + '/kmer_merge345mer.csv'

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
        if kmer == '4mer':
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


def get_ContrastiveLearningDataset(data_path, kmer='4mer', n_views=2, kmer_model_path='empty',
                                   device=torch.device('cpu'),
                                   nokmer=False, cov_meannormalize=False, cov_minmaxnormalize=False, cov_standardization=False,
                                   addvars=False,vars_sqrt=False,kmer_l2_normalize=False, kmerMetric_notl2normalize=False):
    # data_path='/home/wangzy/data/binning/STEC_data/output/for_new_method/data_augmentation/'
    if not data_path.endswith('/'):
        data_path = data_path + '/'
    dataset, namelist = get_kmer_coverage(data_path, kmer, n_views, kmer_model_path, device, nokmer, cov_meannormalize, cov_minmaxnormalize, cov_standardization,addvars,vars_sqrt,kmer_l2_normalize,kmerMetric_notl2normalize)
    return dataset, namelist
