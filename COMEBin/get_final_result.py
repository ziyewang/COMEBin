from scripts.unitem_profile import Profile, make_sure_path_exists
import copy
import os

from collections import defaultdict
import pandas as pd

from scripts.unitem_common import read_bins
from scripts.unitem_markers import Markers
from filter_small_bins import filter_small_bins

# for each bin
# update for checkm marker
def get_binstats(bin_contig_names, markers):
    _, comp, cont = markers.bin_quality(bin_contig_names)

    return comp, cont

def read_bins_nosequences(bin_dirs):
    """Read sequences in bins."""

    bins = defaultdict(lambda: defaultdict(set))

    contigs_in_bins = defaultdict(lambda: {})
    for method_id, bin_dir in bin_dirs.items():
        Header = pd.read_csv(bin_dir, sep='\t', nrows=1)
        cluster_ids = pd.read_csv(bin_dir, sep='\t',header=None, usecols=range(1, Header.shape[1])).values[:, 0]
        namelist = pd.read_csv(bin_dir, sep='\t', header=None, usecols=range(1)).values[:, 0]

        for i in range(len(namelist)):
            bins[method_id][cluster_ids[i]].add(namelist[i])
            contigs_in_bins[namelist[i]][method_id] = cluster_ids[i]

    return bins, contigs_in_bins


# update for checkm marker
def get_bin_quality(orig_bins, methods_sorted, markers):
    # calculate quality of each bin of the orig_bins.
    bin_quality_dict = defaultdict(lambda: {})
    sum_list = []
    sumcont5_list = []


    for method_id in methods_sorted:
        num_5010 = 0
        num_7010 = 0
        num_9010 = 0
        num_505 = 0
        num_705 = 0
        num_905 = 0

        for bin_id in orig_bins[method_id]:
            comp, cont = get_binstats(orig_bins[method_id][bin_id], markers)
            if comp > 50 and cont < 10:
                num_5010 += 1
            if comp > 70 and cont < 10:
                num_7010 += 1
            if comp > 90 and cont < 10:
                num_9010 += 1
            if comp > 50 and cont < 5:
                num_505 += 1
            if comp > 70 and cont < 5:
                num_705 += 1
            if comp > 90 and cont < 5:
                num_905 += 1

        bin_quality_dict[method_id]['num_5010'] = num_5010
        bin_quality_dict[method_id]['num_7010'] = num_7010
        bin_quality_dict[method_id]['num_9010'] = num_9010
        bin_quality_dict[method_id]['num_505'] = num_505
        bin_quality_dict[method_id]['num_705'] = num_705
        bin_quality_dict[method_id]['num_905'] = num_905
        bin_quality_dict[method_id]['sum'] = num_5010 + num_7010 + num_9010 + num_505 + num_705 + num_905
        bin_quality_dict[method_id]['sum_cont5'] = num_505 + num_705 + num_905

        sum_list.append(bin_quality_dict[method_id]['sum'] )
        sumcont5_list.append(bin_quality_dict[method_id]['sum_cont5'])

    sum_max = max(sum_list)
    sum_max_method = []

    sumcont5_remain = []

    for i in range(len(methods_sorted)):
        if sum_list[i] == sum_max:
            sum_max_method.append(methods_sorted[i])
            sumcont5_remain.append(sumcont5_list[i])
    if len(sum_max_method) == 1:
        best_method = sum_max_method[0]
    else:
        best_method = sum_max_method[sumcont5_remain.index(max(sumcont5_remain))]

    return bin_quality_dict, best_method


# update for checkm marker
def savecontigs_with_high_bin_quality(orig_bins, best_method, markers, outpath):
    bin_count_5010 = 0
    bin_count_5005 = 0
    with open(outpath+'/'+best_method+'5010_res.txt','w') as f1:
        with open(outpath+'/'+best_method+'5005_res.txt','w') as f2:
            for bin_id in orig_bins[best_method]:
                comp, cont = get_binstats(orig_bins[best_method][bin_id], markers)
                if comp > 50 and cont < 10:
                    for key in orig_bins[best_method][bin_id]:
                        f1.write(key+'\t'+str(bin_count_5010)+'\n')
                    bin_count_5010 += 1

                if comp > 50 and cont < 5:
                    for key in orig_bins[best_method][bin_id]:
                        f2.write(key+'\t'+str(bin_count_5005)+'\n')
                    bin_count_5005 += 1


def write_estimated_bin_quality(bin_quality_dict, output_file):
    fout = open(os.path.join(output_file), 'w')
    fout.write('Binning_method\tnum_5010\tnum_7010\tnum_9010\tnum_505\tnum_705\tnum_905\tsum\tsum_cont5\n')
    for method_id in bin_quality_dict:
        fout.write(method_id + '\t' + str(bin_quality_dict[method_id]['num_5010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_7010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_9010']) + '\t'
                   + str(bin_quality_dict[method_id]['num_505']) + '\t'
                   + str(bin_quality_dict[method_id]['num_705']) + '\t'
                   + str(bin_quality_dict[method_id]['num_905']) + '\t'
                   + str(bin_quality_dict[method_id]['sum']) + '\t'
                   + str(bin_quality_dict[method_id]['sum_cont5']) + '\n')

    fout.close()


def estimate_bins_quality(bac_mg_table, ar_mg_table, res_path, ignore_kmeans_res = False):
    markers = Markers()

    # bin_dirs = get_bin_dirs(bin_dirs_file)
    filenames = os.listdir(res_path)
    namelist = []
    for filename in filenames:
        if filename.endswith('.tsv'):
            if ignore_kmeans_res:
                if not filename.startswith('weight'):
                    namelist.append(filename)
            else:
                namelist.append(filename)

    namelist.sort()

    bin_dirs = {}
    for res in namelist:
        bin_dirs[res] = (res_path + res + '_bins', 'fa')

    bins, contigs, contigs_in_bins = read_bins(bin_dirs)

    methods_sorted = sorted(bins.keys())
    contig_lens = {cid: len(contigs[cid]) for cid in contigs}
    orig_bins = copy.deepcopy(bins)

    gene_tables = markers.marker_gene_tables(bac_mg_table, ar_mg_table)

    bin_quality_dict, best_method = get_bin_quality(orig_bins, methods_sorted, markers)

    savecontigs_with_high_bin_quality(orig_bins, best_method, markers, res_path)

    output_file = res_path + 'estimate_res.txt'
    write_estimated_bin_quality(bin_quality_dict, output_file)
    return best_method


def estimate_bins_quality_nobins(bac_mg_table, ar_mg_table, res_path, ignore_kmeans_res = False):
    markers = Markers()

    # bin_dirs = get_bin_dirs(bin_dirs_file)
    filenames = os.listdir(res_path)
    namelist = []
    for filename in filenames:
        if filename.endswith('.tsv'):
            if ignore_kmeans_res:
                if not filename.startswith('weight'):
                    namelist.append(filename)
            else:
                namelist.append(filename)

    namelist.sort()

    bin_dirs = {}
    for res in namelist:
        bin_dirs[res] =  (res_path + res)

    bins, contigs = read_bins_nosequences(bin_dirs)

    methods_sorted = sorted(bins.keys())
    contig_lens = {cid: len(contigs[cid]) for cid in contigs}
    orig_bins = copy.deepcopy(bins)

    gene_tables = markers.marker_gene_tables(bac_mg_table, ar_mg_table)

    bin_quality_dict, best_method = get_bin_quality(orig_bins, methods_sorted, markers)

    savecontigs_with_high_bin_quality(orig_bins, best_method, markers, res_path)

    output_file = res_path + 'estimate_res.txt'
    write_estimated_bin_quality(bin_quality_dict, output_file)
    return best_method


def run_get_final_result(logger, args, seed_num, num_threads=40,res_name=None,ignore_kmeans_res=False):
    logger.info("Seed_num:\t" + str(seed_num))

    if not (args.bac_mg_table and args.ar_mg_table):
        logger.info("Run unitem profile:\t" + str(seed_num))
        # run unitem_profile
        bin_dirs = {}
        if res_name==None:
            res_name = 'weight_partialseed_kmeans_algofull_rand_1_k_' + str(seed_num + 1) + '_result.tsv'
        bin_dirs[res_name] = (args.output_path + '/cluster_res/' + res_name + '_bins', 'fa')

        output_dir = args.output_path + '/cluster_res/unitem_profile'

        if not (os.path.exists(output_dir)):
            make_sure_path_exists(output_dir)
            profile = Profile(num_threads)
            profile.run(bin_dirs,
                        output_dir)

        bac_mg_table = output_dir + '/binning_methods/' + res_name + '/checkm_bac/marker_gene_table.tsv'
        ar_mg_table = output_dir + '/binning_methods/' + res_name + '/checkm_ar/marker_gene_table.tsv'
    else:
        bac_mg_table = args.bac_mg_table
        ar_mg_table = args.ar_mg_table

    # best_method = estimate_bins_quality(bac_mg_table, ar_mg_table, args.output_path + '/cluster_res/',ignore_kmeans_res=ignore_kmeans_res)
    best_method = estimate_bins_quality_nobins(bac_mg_table, ar_mg_table, args.output_path + '/cluster_res/',ignore_kmeans_res=ignore_kmeans_res)

    logger.info('Final result:\t'+args.output_path + '/cluster_res/'+best_method)
    filter_small_bins(args.contig_file, args.output_path + '/cluster_res/'+best_method, args)

