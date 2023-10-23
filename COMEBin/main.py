import argparse
import logging
import os
import pandas as pd

from comebin_version import __version__ as ver
from train_CLmodel import train_CLmodel
from cluster import cluster


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
                        help='path to dataset')
    CLtraining_subparsers.add_argument('--kmer_model_path', metavar='DIR', default='empty',
                        help='kmer_model_path')
    CLtraining_subparsers.add_argument('--output_path', metavar='DIR', default='output',
                                       help='Output path.')
    CLtraining_subparsers.add_argument('-j', '--workers', default=5, type=int, metavar='N',
                        help='number of data loading workers (default: 5)')
    CLtraining_subparsers.add_argument('--epochs', default=200, type=int, metavar='N',
                        help='number of total epochs to run')
    CLtraining_subparsers.add_argument('--covmodelepochs', default=50, type=int, metavar='N',
                        help='number of total epochs to pretrain coverage model')

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
    CLtraining_subparsers.add_argument('--add_model_for_coverage', action="store_true",
                        help='add_model_for_coverage.')

    CLtraining_subparsers.add_argument('--nokmer', action="store_true",
                        help='only use coverage information.')

    CLtraining_subparsers.add_argument('--pretrain_coveragemodel', action="store_true",
                        help='train ')

    CLtraining_subparsers.add_argument('--emb_szs_forcov', default=2048, type=int,
                        help='embedding size for hidden layer (default: 2048)')
    CLtraining_subparsers.add_argument('--out_dim_forcov', default=128, type=int,
                        help='embedding size for hidden layer (default: 128)')
    CLtraining_subparsers.add_argument('--n_layer_forcov', default=3, type=int,
                        help='n_layer_forcov (default: 3)')

    CLtraining_subparsers.add_argument('--earlystop', action="store_true",
                                       help='earlystop.')

    CLtraining_subparsers.add_argument('--pretrain_kmer_model_path', metavar='DIR', default='no',
                        help='pretrain_kmer_model_path')

    CLtraining_subparsers.add_argument('--not_load_kmermetric_state', action="store_true",
                        help='not_load_kmermetric_state')

    CLtraining_subparsers.add_argument('--finetunepretrainmodel', action="store_true",
                        help='finetunepretrainmodel')
    CLtraining_subparsers.add_argument('--finetunelr_ratio', default=0.1, type=float, help='finetune learning rate')

    CLtraining_subparsers.add_argument('--addvars', action="store_true",
                                       help='addvars')
    CLtraining_subparsers.add_argument('--vars_sqrt', action="store_true",
                                       help='vars_sqrt')

    CLtraining_subparsers.add_argument('--log-every-n-steps', default=20, type=int,
                        help='Log every n steps')

    CLtraining_subparsers.add_argument('--temperature', default=0.1, type=float,
                        help='softmax temperature (default: 0.1)')
    CLtraining_subparsers.add_argument('--covmodel_temperature', default=0.1, type=float,
                        help='softmax temperature (default: 0.1)')

    CLtraining_subparsers.add_argument('--n_views', default=6, type=int, metavar='N',
                        help='Number of views for contrastive learning training.')
    CLtraining_subparsers.add_argument('--contig_len', default = 1000, type=int, metavar='N',
                        help='mininum contig length for training')
    CLtraining_subparsers.add_argument('--notuse_scheduler', action='store_true',
                        help='notuse_scheduler')
    CLtraining_subparsers.add_argument('--fp16-precision', action='store_true',
                        help='Whether or not to use 16-bit precision GPU training.')
    CLtraining_subparsers.add_argument('--addcovloss', action='store_true',
                        help='addcovloss')
    CLtraining_subparsers.add_argument('--lambdaloss2', default=1, type=float,
                        help='lambdaloss2 (default: 1)')

    CLtraining_subparsers.add_argument('--addkmerloss', action='store_true',
                                       help='addkmerloss')
    CLtraining_subparsers.add_argument('--lambdakmerloss2', default=1, type=float,
                                       help='lambdakmerloss2 (default: 1)')
    CLtraining_subparsers.add_argument('--kmermodel_temperature', default=0.1, type=float,
                                       help='kmermodel_temperature softmax temperature (default: 0.1)')

    CLtraining_subparsers.add_argument('--cov_meannormalize', action='store_true',
                                       help='cov_meannormalize')
    CLtraining_subparsers.add_argument('--cov_minmaxnormalize', action='store_true',
                                       help='cov_minmaxnormalize')
    CLtraining_subparsers.add_argument('--cov_standardization', action='store_true',
                                       help='cov_standardization')
    CLtraining_subparsers.add_argument('--notuseaug0', action='store_true',
                                       help='notuseaug0')
    CLtraining_subparsers.add_argument('--lossbalance', action='store_true',
                                       help='lossbalance ((whole_loss[epoch0]/covloss[epoch0])')

    CLtraining_subparsers.add_argument('--kmer_l2_normalize', action='store_true',
                                       help='kmer_l2_normalize (used for nokmerMetric)')

    CLtraining_subparsers.add_argument('--kmerMetric_notl2normalize', action='store_true',
                                       help='kmerMetric_notl2normalize (used for kmerMetric)')
    CLtraining_subparsers.add_argument('--covmodel_notl2normalize', action='store_true',
                                       help='covmodel_notl2normalize (used for covmodel)')




    #############################################################################################
    ############################################ cluster NoContrast #####################################
    ### Command-line arguments and options for running the COMEBin using the original features.

    NoContrast_subparsers = subparsers.add_parser('nocontrast',
                                                   help='Cluster the contigs using original features.')

    NoContrast_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    NoContrast_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    NoContrast_subparsers.add_argument('--data', metavar='DIR', default='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/',
                                       help='path to dataset')
    NoContrast_subparsers.add_argument('--output_path', type=str, default='temp_output', help=("The output path"))

    NoContrast_subparsers.add_argument('--cluster_num', default=0, type=int,
                                       help='Add cluster number to run partial seed method (default: 0)')

    NoContrast_subparsers.add_argument('--not_run_infomap', action='store_true',
                                       help='Do not run infomap.')

    NoContrast_subparsers.add_argument('--not_l2normaize', action='store_true',
                                       help='Do not run l2normaize for embeddings.')

    NoContrast_subparsers.add_argument('--contig_len', default = 1001, type=int, metavar='N',
                                       help='mininum contig length for clustering')



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
                                              help='num_threads for generating augmentation data.')

    #############################################################################################
    ############################################ cluster #####################################
    ### Command-line arguments and options for running the Leiden-based clustering.

    clustering_subparsers = subparsers.add_parser('bin',
                                                  help='Cluster the contigs.')

    clustering_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    clustering_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    clustering_subparsers.add_argument('--emb_file', type=str, help=("The embedding feature file."))
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


    #############################################################################################
    ############################################ get final results #####################################
    ### Command-line arguments and options for generating the final binning result from the Leiden clustering results.

    get_result_subparsers = subparsers.add_parser('get_result',
                                                  help='Generate the final results from the Leiden clustering results.')

    get_result_subparsers.add_argument('--contig_file', type=str, help=("The contigs file."))
    get_result_subparsers.add_argument('--seed_file', type=str, help=("The marker seed file."))
    get_result_subparsers.add_argument('--emb_file', type=str, help=("The embedding feature file."))
    get_result_subparsers.add_argument('--output_path', type=str, help=("The output path"))
    get_result_subparsers.add_argument('--binning_res_path', type=str, help=("The path to get Leiden clustering results"))
    get_result_subparsers.add_argument('--contig_len', default = 1001, type=int, metavar='N',
                                       help='mininum contig length for clustering')
    get_result_subparsers.add_argument('--num_threads', default=10, type=int,
                                       help='num_threads for getting final result.')

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

    os.makedirs(args.output_path, exist_ok=True)
    handler = logging.FileHandler(args.output_path+'/comebin.log')
    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    ## training
    if args.subcmd == 'train':
        logger.info('train')
        train_CLmodel(logger,args)

    ## clustering
    if args.subcmd == 'bin':
        logger.info('bin')
        from utils import gen_seed

        num_threads = args.num_threads
        _ = gen_seed(logger, args.contig_file, num_threads, args.contig_len, marker_name="bacar_marker", quarter="2quarter")

        cluster(logger, args)


    ## clustering NoContrast
    if args.subcmd == 'nocontrast':
        logger.info('NoContrast mode')
        from utils import get_kmer_coverage_aug0

        X_t, covMat, compositMat, namelist = get_kmer_coverage_aug0(args.data)

        X_t_df = pd.DataFrame(X_t, index=namelist)
        os.makedirs(args.output_path+'/combine_novars', exist_ok=True)
        outfile = args.output_path+'/combine_novars/combine_feature.tsv'
        X_t_df.to_csv(outfile, sep='\t', header=True)

        covMat_df = pd.DataFrame(covMat, index=namelist)
        os.makedirs(args.output_path+'/covMat', exist_ok=True)
        outfile = args.output_path+'/covMat/covMat_feature.tsv'
        covMat_df.to_csv(outfile, sep='\t', header=True)

        compositMat_df = pd.DataFrame(compositMat, index=namelist)
        os.makedirs(args.output_path+'/compositMat', exist_ok=True)
        outfile = args.output_path+'/compositMat/compositMat_feature.tsv'
        compositMat_df.to_csv(outfile, sep='\t', header=True)

        logger.info('NoContrast mode: generate features (aug0)')
        ori_outpath = args.output_path

        logger.info('NoContrast mode (combine) bin')
        args.output_path = ori_outpath +'/combine_novars'
        args.emb_file = args.output_path+'/combine_feature.tsv'
        cluster(logger,args)

        logger.info('NoContrast mode (coverage) bin')
        args.output_path = ori_outpath +'/covMat'
        args.emb_file = args.output_path+'/covMat_feature.tsv'
        cluster(logger,args)

        logger.info('NoContrast mode (kmer) bin')
        args.output_path = ori_outpath +'/compositMat'
        args.emb_file = args.output_path+'/compositMat_feature.tsv'
        cluster(logger,args)



    ##### generate_aug_data fastafile
    if args.subcmd == 'generate_aug_data':
        logger.info('generate_aug_data: fastafile')

        from data_aug.generate_augfasta_and_saveindex import run_gen_augfasta
        from data_aug.gen_cov import run_gen_cov
        from data_aug.gen_var import run_gen_cov_var

        run_gen_augfasta(logger, args)
        run_gen_cov(logger, args)
        run_gen_cov_var(logger, args)

    ###Generate the final results from the Leiden clustering results
    if args.subcmd == 'get_result':
        logger.info('get_result')
        from utils import gen_seed
        from get_final_result import run_get_final_result

        num_threads = args.num_threads
        seed_num = gen_seed(logger, args.contig_file, num_threads, args.contig_len, marker_name="bacar_marker", quarter="2quarter")

        run_get_final_result(logger, args, seed_num, num_threads, ignore_kmeans_res=True)


if __name__ == '__main__':
    main()

