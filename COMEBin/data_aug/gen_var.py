# modified from https://github.com/BigDataBiology/SemiBin/blob/main/SemiBin/generate_coverage.py
import multiprocessing
import traceback
from multiprocessing.pool import Pool

import os
from atomicwrites import atomic_write
import pandas as pd
import numpy as np
from itertools import groupby
from typing import List, Optional, Union, Dict, Tuple

### Return error message when using multiprocessing
def error(msg, *args):
    return multiprocessing.get_logger().error(msg, *args)


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)
        except Exception as e:
            error(traceback.format_exc())
            raise

        return result


class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)


def _checkback(msg):
    msg[1].info('Processed:{}'.format(msg[0]))


# modifed from https://github.com/BigDataBiology/SemiBin/blob/3bad22c58e710d8a5455f7411bc8d4202d557c61/SemiBin/generate_coverage.py#L5
def calculate_coverage_var_samplebyindex(depth_file: str, augpredix: str, aug_seq_info_dict: Dict[str, Tuple[int, int]],
                                         logger, edge: int = 0, contig_threshold: int = 1000):
    """
    Calculate coverage variance per contig and save the results to a CSV file.

    :param depth_file: Path to the position depth file generated from bedtools genomecov.
    :param augpredix: A prefix used in the output file name.
    :param aug_seq_info_dict: A dictionary containing contig information (start and end positions) as (start, end) tuples.
    :param edge: The number of bases to exclude from the edges of each contig (default is 0).
    :param contig_threshold: The minimum depth threshold for a contig to be considered (default is 1000).

    :return: A tuple containing the path to the processed depth file and the logger object.
    """

    contigs = []
    var_coverage = []

    for contig_name, lines in groupby(open(depth_file), lambda ell: ell.split('\t', 1)[0]):
        depth_value = []
        for line in lines:
            line_split = line.strip().split('\t')
            length = int(float(line_split[2])) - int(float(line_split[1]))
            value = int(float(line_split[3]))
            depth_value.extend([value] * length)

        cov_threshold = contig_threshold

        if len(depth_value) <= cov_threshold:
            continue
        start = aug_seq_info_dict[contig_name][0]
        end = aug_seq_info_dict[contig_name][1]

        depth_value_ = depth_value[start + edge:end + 1 - edge]

        var = np.var(depth_value_)
        var_coverage.append(var)
        contigs.append(contig_name)

    contig_cov = pd.DataFrame(
        {'{0}_var'.format(depth_file): var_coverage,
         }, index=contigs)

    with atomic_write(depth_file + '_' + augpredix + '_data_var.csv', overwrite=True) as ofile:
        contig_cov.to_csv(ofile, sep='\t')

    return (depth_file, logger)


def calculate_coverage_var(depth_file: str, logger, edge: int = 0, contig_threshold: int = 1000, sep: Optional[str] = None,
                           contig_threshold_dict: Optional[Dict[str, int]] = None):
    """
    Calculate coverage variance per contig and save the results to a CSV file.

    :param depth_file: Path to the position depth file generated from bedtools genomecov.
    :param edge: The number of bases to exclude from the edges of each contig (default is 0).
    :param contig_threshold: The minimum depth threshold for a contig to be considered (default is 1000).
    :param sep: Separator for distinguishing sample names in contig names (default is None).
    :param contig_threshold_dict: A dictionary containing sample-specific contig thresholds when `sep` is provided (default is None).

    :return: A tuple containing the path to the processed depth file and the logger object.
    """
    contigs = []
    var_coverage = []

    for contig_name, lines in groupby(open(depth_file), lambda ell: ell.split('\t', 1)[0]):
        depth_value = []
        for line in lines:
            line_split = line.strip().split('\t')
            length = int(float(line_split[2])) - int(float(line_split[1]))
            value = int(float(line_split[3]))
            depth_value.extend([value] * length)

        if sep is None:
            cov_threshold = contig_threshold
        else:
            sample_name = contig_name.split(sep)[0]
            cov_threshold = contig_threshold_dict[sample_name]
        if len(depth_value) <= cov_threshold:
            continue
        # depth_value_ = depth_value[edge:-edge]
        depth_value_ = depth_value

        # modified
        var = np.var(depth_value_)
        var_coverage.append(var)
        contigs.append(contig_name)

    contig_cov = pd.DataFrame(
        {'{0}_var'.format(depth_file): var_coverage,
         }, index=contigs)

    with atomic_write(depth_file + '_aug0_data_var.csv', overwrite=True) as ofile:
        contig_cov.to_csv(ofile, sep='\t')

    return (depth_file, logger)


def gen_cov_var_from_bedout(logger, out_path, depth_file_path, num_process=10, num_aug=5,edge=0, contig_len=1000):
    filenames = os.listdir(depth_file_path)
    namelist = []
    for filename in filenames:
        if filename.endswith('_depth.txt'):
            namelist.append(filename)

    namelist.sort()

    pool = LoggingPool(num_process) if num_process != 0 else LoggingPool()
    ##generate coverage for original data
    for i in range(len(namelist)):
        depth_file = depth_file_path + namelist[i]
        pool.apply_async(
            calculate_coverage_var,
            args=(depth_file, logger, edge, contig_len),
            callback=_checkback)

    pool.close()
    pool.join()

    # merge coverage files
    for nameid in range(len(namelist)):
        # cov_file = depth_file_path + namelist[
        #     nameid] + '_aug0_data_cov_edge75.csv'
        cov_file = depth_file_path + namelist[
            nameid] + '_aug0_data_var.csv'
        res_mat = pd.read_csv(cov_file, sep='\t', header=0, index_col=0)
        if nameid == 0:
            joined = res_mat
        else:
            joined = res_mat.join(joined, how="inner")

    # outfile = out_path + 'aug0_datacoverage_mean_edge75.tsv'
    outfile = out_path + 'aug0_datacoverage_var.tsv'
    joined.to_csv(outfile, sep='\t', header=True)

    for i in range(num_aug):
        outdir = out_path + 'aug' + str(i + 1)
        aug_seq_info_out_file = outdir + '/sequences_aug' + str(i + 1) + '.fasta' + '.aug_seq_info.tsv'
        aug_seq_info_dict = read_aug_seq_info(aug_seq_info_out_file)

        # num_process = 40
        pool = LoggingPool(num_process) if num_process != 0 else LoggingPool()

        ####generate coverage files
        for nameid in range(len(namelist)):
            depth_file = depth_file_path + namelist[nameid]
            pool.apply_async(
                calculate_coverage_var_samplebyindex,
                args=(depth_file, 'aug' + str(i + 1), aug_seq_info_dict, logger, edge, contig_len),
                callback=_checkback)

        pool.close()
        pool.join()

        # merge coverage files
        for nameid in range(len(namelist)):
            # cov_file = depth_file_path + namelist[
            #     nameid] + '_aug' + str(i + 1) + '_data_cov_edge75.csv'
            cov_file = depth_file_path + namelist[
                nameid] + '_aug' + str(i + 1) + '_data_var.csv'

            res_mat = pd.read_csv(cov_file, sep='\t', header=0, index_col=0)
            if nameid == 0:
                joined = res_mat
            else:
                joined = res_mat.join(joined, how="inner")

        logger.info("Finish calculating coverage variance for aug_"+str(i+1))
        outfile = out_path + 'aug' + str(i + 1) + '_datacoverage_var.tsv'
        joined.to_csv(outfile, sep='\t', header=True)


def read_aug_seq_info(aug_seq_info_out_file):
    aug_seq_info = pd.read_csv(aug_seq_info_out_file, sep='\t', header=0).values[:]
    aug_seq_info_dict = {}
    for i in range(len(aug_seq_info)):
        aug_seq_info_dict[aug_seq_info[i][0]] = [aug_seq_info[i][1], aug_seq_info[i][2]]

    return aug_seq_info_dict


def run_gen_cov_var(logger, args):
    logger.info("Generate coverage variance files from bam files.")

    if not args.out_augdata_path.endswith('/'):
        args.out_augdata_path = args.out_augdata_path + '/'

    out = args.out_augdata_path + 'depth/'
    gen_cov_var_from_bedout(logger, args.out_augdata_path, out, num_aug=args.n_views-1,contig_len=args.contig_len,num_process=args.num_threads)
