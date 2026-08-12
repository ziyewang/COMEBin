import argparse
import copy
import faulthandler
import logging
import os
from pathlib import Path
import numpy as np

from comebin_version import __version__ as ver
from train_CLmodel import train_CLmodel
from cluster import cluster
from utils import DEFAULT_HMM_EVALUE, DEFAULT_RANDOM_SEED, configure_reproducibility


def arguments():
    """
    COMEBin: A contig binning method based on Contrastive Multi-view representation learning.

    :return: Parsed command-line arguments
    """
    doc = f"""COMEBin: a contig binning method based on COntrastive Multi-viEw representation learning.
    Version:{ver}"""
    parser = argparse.ArgumentParser(
        prog="comebin",
        description=doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage="%(prog)s outdir contig_file path_to_bamfiles [options]")
    parser.version = ver

    parser.add_argument('-v','--version',
                        action='version',
                        help='COMEBin version')

    subparsers = parser.add_subparsers(title='COMEBin subcommands',
                                       dest='subcmd',
                                       metavar='')

    #############################################################################################
    ############################################ CLtraining #####################################
    ### Command-line arguments and options for training the network.

    CLtraining_subparsers = subparsers.add_parser('train',
                                                  help='Train the model based on the augmentation data')

    CLtraining_subparsers.add_argument('--data', metavar='DIR', default='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/',
                        help='path to a completed binary augmentation and coverage bundle')
    CLtraining_subparsers.add_argument('--output_path', metavar='DIR', default='output',
                                       help='Output path.')
    CLtraining_subparsers.add_argument('--epochs', default=200, type=int, metavar='N',
                        help='number of total epochs to run')

    CLtraining_subparsers.add_argument('-b', '--batch_size', default=1024, type=int,
                        metavar='N',
                        help='mini-batch size (default: 1024), this is the total '
                             'batch size of all GPUs on the current node when '
                             'using Data Parallel or Distributed Data Parallel')
    CLtraining_subparsers.add_argument('--lr', '--learning-rate', default=0.001, type=float,
                        metavar='LR', help='initial learning rate', dest='lr')
    CLtraining_subparsers.add_argument('--wd', '--weight-decay', default=1e-4, type=float,
                        metavar='W', help='weight decay (default: 1e-4)',
                        dest='weight_decay')
    CLtraining_subparsers.add_argument('--out_dim', default=128, type=int,
                        help='feature dimension (default: 128)')

    CLtraining_subparsers.add_argument('--dropout_value', default=0.2, type=float,
                        help='dropout_value (default: 0.2)')
    CLtraining_subparsers.add_argument('--emb_szs', default=2048, type=int,
                        help='embedding size for hidden layer (default: 2048)')
    CLtraining_subparsers.add_argument('--n_layer', default=3, type=int,
                        help='n layers (default: 3)')

    CLtraining_subparsers.add_argument('--emb_szs_forcov', default=2048, type=int,
                        help='embedding size for hidden layer (default: 2048)')
    CLtraining_subparsers.add_argument('--out_dim_forcov', default=128, type=int,
                        help='embedding size for hidden layer (default: 128)')
    CLtraining_subparsers.add_argument('--n_layer_forcov', default=3, type=int,
                        help='n_layer_forcov (default: 3)')

    CLtraining_subparsers.add_argument('--earlystop', action="store_true",
                                       help='earlystop.')

    CLtraining_subparsers.add_argument('--addvars', action="store_true",
                                       help='addvars')
    CLtraining_subparsers.add_argument('--vars_sqrt', action="store_true",
                                       help='vars_sqrt')

    CLtraining_subparsers.add_argument('--log-every-n-steps', default=20, type=int,
                        help='Log every n steps')

    CLtraining_subparsers.add_argument(
        '--temperature', default=None, type=float,
        help='InfoNCE temperature; when omitted, use 0.07 for assembly N50 > 10000 and 0.15 otherwise.')

    CLtraining_subparsers.add_argument('--n_views', default=6, type=int, metavar='N',
                        help='Number of views for contrastive learning training.')
    CLtraining_subparsers.add_argument('--contig_len', default = 1000, type=int, metavar='N',
                        help='mininum contig length for training')
    CLtraining_subparsers.add_argument('--notuse_scheduler', action='store_true',
                        help='notuse_scheduler')
    CLtraining_subparsers.add_argument('--fp16-precision', action='store_true',
                        help='Whether or not to use 16-bit precision GPU training.')
    CLtraining_subparsers.add_argument('--device', default='cuda',
                                       help='NN training device: cuda, cuda:N, or cpu.')

    CLtraining_subparsers.add_argument('--cov_meannormalize', action='store_true',
                                       help='cov_meannormalize')
    CLtraining_subparsers.add_argument('--cov_minmaxnormalize', action='store_true',
                                       help='cov_minmaxnormalize')
    CLtraining_subparsers.add_argument('--cov_standardization', action='store_true',
                                       help='cov_standardization')
    CLtraining_subparsers.add_argument('--notuseaug0', action='store_true',
                                       help='notuseaug0')
    CLtraining_subparsers.add_argument('--num_threads', default=10, type=int,
                                       help='PyTorch intra-op threads; does not control DataLoader or native BLAS/OpenMP pools.')
    CLtraining_subparsers.add_argument('--seed', default=DEFAULT_RANDOM_SEED, type=int,
                                       help='optional random seed for reproducible training.')



    #############################################################################################
    ############################################ cluster NoContrast #####################################
    ### Command-line arguments and options for running the COMEBin using the original features.

    NoContrast_subparsers = subparsers.add_parser('nocontrast',
                                                   help='Cluster the contigs using original features.')

    NoContrast_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    NoContrast_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    NoContrast_subparsers.add_argument('--data', metavar='DIR', default='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/',
                                       help='path to a completed binary augmentation and coverage bundle')
    NoContrast_subparsers.add_argument('--output_path', type=str, default='temp_output', help=("The output path"))

    NoContrast_subparsers.add_argument('--cluster_num', default=0, type=int,
                                       help='Add cluster number to run partial seed method (default: 0)')

    NoContrast_subparsers.add_argument('--not_run_infomap', action='store_true',
                                       help='Do not run infomap.')

    NoContrast_subparsers.add_argument('--not_l2normaize', action='store_true',
                                       help='Do not run l2normaize for embeddings.')

    NoContrast_subparsers.add_argument('--contig_len', default = 1001, type=int, metavar='N',
                                       help='mininum contig length for clustering')
    NoContrast_subparsers.add_argument('--num_threads', default=10, type=int,
                                       help='num_threads for binning.')
    NoContrast_subparsers.add_argument('--leiden_workers', default=None, type=int,
                                       help='concurrent Leiden worker processes; defaults to num_threads.')
    NoContrast_subparsers.add_argument('--max_edges', default=100, type=int,
                                       help='HNSW neighbors retained per contig for Leiden graph construction.')
    NoContrast_subparsers.add_argument('--seed', default=DEFAULT_RANDOM_SEED, type=int,
                                       help='optional random seed for reproducible clustering.')



    #############################################################################################
    ############################################ generate aug data #####################################
    ### Command-line arguments and options for data augmentation.

    generate_aug_data_subparsers = subparsers.add_parser('generate_aug_data',
                                                      help='Generate the augmentation data and features from the fasta file and bam files.')
    generate_aug_data_subparsers.add_argument('--contig_file', type=str, help=("The original contigs file."))
    generate_aug_data_subparsers.add_argument('--out_augdata_path', type=str, help=("The output path to save the augmentation data"))
    generate_aug_data_subparsers.add_argument('--n_views', default=6, type=int,
                                           help='n_views for generating augmentation data.')
    generate_aug_data_subparsers.add_argument('--bam_file_path', type=str, help=("The path to access the bam files."))


    generate_aug_data_subparsers.add_argument('--contig_len', default = 1000, type=int, metavar='N',
                                       help='mininum contig length for augmentation')
    generate_aug_data_subparsers.add_argument('--num_threads', default=10, type=int,
                                              help='concurrent coverage workers, capped by the number of BAM files.')
    generate_aug_data_subparsers.add_argument('--seed', default=DEFAULT_RANDOM_SEED, type=int,
                                              help='optional random seed for reproducible data augmentation.')

    #############################################################################################
    ############################################ cluster #####################################
    ### Command-line arguments and options for running the Leiden-based clustering.

    clustering_subparsers = subparsers.add_parser('bin',
                                                  help='Cluster the contigs.')

    clustering_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    clustering_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    clustering_subparsers.add_argument('--emb_file', type=str,
                                       help=("Canonical embeddings.npy."))
    clustering_subparsers.add_argument('--output_path', type=str, help=("The output path"))

    clustering_subparsers.add_argument('--cluster_num', default=0, type=int,
                                       help='Add cluster number to run partial seed method (default: 0)')

    clustering_subparsers.add_argument('--not_run_infomap', action='store_true',
                                       help='Do not run infomap.')
    clustering_subparsers.add_argument('--not_l2normaize', action='store_true',
                                       help='Do not run l2normaize for embeddings.')


    clustering_subparsers.add_argument('--contig_len', default = 1001, type=int, metavar='N',
                                       help='mininum contig length for clustering')
    clustering_subparsers.add_argument('--num_threads', default=10, type=int,
                                              help='num_threads for binning.')
    clustering_subparsers.add_argument('--leiden_workers', default=None, type=int,
                                       help='concurrent Leiden worker processes; defaults to num_threads.')
    clustering_subparsers.add_argument('--max_edges', default=100, type=int,
                                       help='HNSW neighbors retained per contig for Leiden graph construction.')
    clustering_subparsers.add_argument('--hmm_evalue', default=DEFAULT_HMM_EVALUE, type=float,
                                       help='HMMER sequence E-value used when marker HMMs do not define TC cutoffs.')
    clustering_subparsers.add_argument('--seed', default=DEFAULT_RANDOM_SEED, type=int,
                                       help='optional random seed for reproducible clustering.')


    #############################################################################################
    ############################################ get final results #####################################
    ### Command-line arguments and options for generating the final binning result from the Leiden clustering results.

    get_result_subparsers = subparsers.add_parser('get_result',
                                                  help='Generate the final results from the Leiden clustering results.')

    get_result_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    get_result_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    get_result_subparsers.add_argument('--emb_file', type=str,
                                       help=("Canonical embeddings.npy used to resolve the shared result axis."))
    get_result_subparsers.add_argument('--output_path', type=str, help=("The output path"))
    get_result_subparsers.add_argument('--binning_res_path', type=str, help=("The path to get Leiden clustering results"))
    get_result_subparsers.add_argument('--contig_len', default = 1001, type=int, metavar='N',
                                       help='mininum contig length for clustering')
    get_result_subparsers.add_argument('--num_threads', default=10, type=int,
                                       help='num_threads for getting final result.')
    get_result_subparsers.add_argument('--max_edges', default=100, type=int,
                                       help='HNSW neighbor count used by the completed Leiden run.')
    get_result_subparsers.add_argument('--hmm_evalue', default=DEFAULT_HMM_EVALUE, type=float,
                                       help='HMMER sequence E-value used when marker HMMs do not define TC cutoffs.')
    get_result_subparsers.add_argument('--seed', default=DEFAULT_RANDOM_SEED, type=int,
                                       help='optional random seed for reproducible runs.')

    get_result_subparsers.add_argument('--bac_mg_table', type=str, help=("bac_mg_table (bacteria marker gene information)"))
    get_result_subparsers.add_argument('--ar_mg_table', type=str, help=("ar_mg_table (archea marker gene information)"))

    args = parser.parse_args()
    return args



def main():
    """
    The main function of the COMEBin program.

    Functionality:
        - Initializes logging for the program.
        - Executes different subcommands based on user input.
        - Subcommands include: 'train', 'bin', 'nocontrast', 'generate_aug_data', and 'get_result'.
        - Subcommands perform various tasks such as data augmentation, training and clustering.
    """
    args = arguments()

    # logging
    logger = logging.getLogger('COMEBin\t'+ver)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(message)s')
    console_hdr = logging.StreamHandler()
    console_hdr.setFormatter(formatter)

    logger.addHandler(console_hdr)

    if args.subcmd == 'generate_aug_data':
        args.output_path = args.out_augdata_path

    output_path = Path(args.output_path).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    args.output_path = str(output_path)
    handler = logging.FileHandler(output_path / 'comebin.log')
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    faulthandler.enable(file=handler.stream, all_threads=True)
    if args.seed is not None:
        logger.info("Configure reproducibility: start with seed=%d, pid=%d.", args.seed, os.getpid())
        configure_reproducibility(args.seed)
        logger.info("Configure reproducibility: finished; seed=%d, PYTHONHASHSEED=%s, CUBLAS_WORKSPACE_CONFIG=%s.",
                    args.seed, os.environ.get('PYTHONHASHSEED', 'not set'),
                    os.environ.get('CUBLAS_WORKSPACE_CONFIG', 'not set'))

    ## training
    if args.subcmd == 'train':
        logger.info('train')
        try:
            train_CLmodel(logger, args)
        except Exception:
            logger.exception("Training failed.")
            raise

    ## clustering
    if args.subcmd == 'bin':
        logger.info('bin')
        from utils import gen_seed

        num_threads = args.num_threads
        _ = gen_seed(logger, args.contig_file, num_threads, args.contig_len, marker_name="bacar_marker",
                     quarter="2quarter", hmm_evalue=args.hmm_evalue)

        cluster(logger, args)


    ## clustering NoContrast
    if args.subcmd == 'nocontrast':
        logger.info('NoContrast mode')
        from cluster import cluster_features, prepare_cluster_axis
        from get_augfeature import load_nocontrast_view0

        coverage, kmer, contig_ids, contig_lengths = load_nocontrast_view0(
            args.data, logger)
        retained, contig_ids, contig_lengths, seed_indexes, seed_count, n50 = \
            prepare_cluster_axis(
                contig_ids, contig_lengths, args.seed_file,
                args.contig_len, logger)
        logger.info(
            'Prepared shared nocontrast axis: retained_contigs=%d, '
            'unique_seeds=%d, retained_seeds=%d.',
            len(contig_ids), seed_count, len(seed_indexes))
        if not retained.all():
            coverage = coverage[retained]
            kmer = kmer[retained]

        original_output = Path(args.output_path)
        combined = np.empty(
            (len(contig_ids), coverage.shape[1] + kmer.shape[1]),
            dtype=np.float32)
        combined[:, :coverage.shape[1]] = coverage
        combined[:, coverage.shape[1]:] = kmer

        logger.info('NoContrast mode (combine) bin')
        combined_args = copy.copy(args)
        combined_args.output_path = str(original_output / 'combine_novars')
        cluster_features(
            logger, combined_args, combined, contig_ids, contig_lengths,
            seed_indexes, seed_count, n50)
        del combined

        logger.info('NoContrast mode (coverage) bin')
        coverage_args = copy.copy(args)
        coverage_args.output_path = str(original_output / 'covMat')
        cluster_features(
            logger, coverage_args, coverage, contig_ids, contig_lengths,
            seed_indexes, seed_count, n50)
        del coverage

        logger.info('NoContrast mode (kmer) bin')
        kmer_args = copy.copy(args)
        kmer_args.output_path = str(original_output / 'compositMat')
        cluster_features(
            logger, kmer_args, kmer, contig_ids, contig_lengths,
            seed_indexes, seed_count, n50)
        del kmer



    ##### generate_aug_data fastafile
    if args.subcmd == 'generate_aug_data':
        logger.info('generate_aug_data: fastafile')

        from data_aug.generate_augfasta_and_saveindex import run_gen_augfasta
        from data_aug.gen_cov import run_gen_coverage

        run_gen_augfasta(logger, args)
        run_gen_coverage(logger, args)

    ###Generate the final results from the Leiden clustering results
    if args.subcmd == 'get_result':
        logger.info('get_result')
        from utils import gen_seed
        from get_final_result import run_get_final_result

        num_threads = args.num_threads
        seed_num = gen_seed(logger, args.contig_file, num_threads, args.contig_len, marker_name="bacar_marker",
                            quarter="2quarter", hmm_evalue=args.hmm_evalue)

        run_get_final_result(logger, args, seed_num, num_threads, ignore_kmeans_res=True)


if __name__ == '__main__':
    main()
