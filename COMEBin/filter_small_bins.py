import gzip
import os
import shutil


def gen_bins(fastafile: str, resultfile: str, outputdir: str) -> None:
    """
    Generate bins from contigs based on a result file and save them to the specified output directory.

    :param fastafile: The path to the input FASTA file containing contigs.
    :param resultfile: The path to the result file that associates contigs with clusters.
    :param outputdir: The output directory where bins will be saved.
    :return: None
    """
    print("Processing file:\t{}".format(fastafile))
    sequences = {}
    if fastafile.endswith("gz"):
        with gzip.open(fastafile, 'r') as f:
            for line in f:
                line = str(line, encoding="utf-8")
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    else:
        with open(fastafile, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if " " in line:
                        seq, others = line.split(' ', 1)
                        sequences[seq] = ""
                    else:
                        seq = line.rstrip("\n")
                        sequences[seq] = ""
                else:
                    sequences[seq] += line.rstrip("\n")
    print("Reading Map:\t{}".format(resultfile))
    dic = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
            except:
                dic[cluster_name] = []
                dic[cluster_name].append(contig_name)
    print("Writing bins:\t{}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    bin_name = 0
    for _, cluster in dic.items():
        binfile = os.path.join(outputdir, "{}.fa".format(bin_name))
        with open(binfile, "w") as f:
            for contig_name in cluster:
                contig_name = ">" + contig_name
                try:
                    sequence = sequences[contig_name]
                except:
                    bin_name += 1
                    continue
                f.write(contig_name + "\n")
                f.write(sequence + "\n")
                bin_name += 1


def filter_small_bins(logger, fastafile: str, resultfile: str, args, minbinsize: int = 200000) -> None:
    """
    Filter small bins from the result file.

    :param fastafile: The path to the input FASTA file containing contigs.
    :param resultfile: The path to the binning result file.
    :param args: The additional arguments used in the process.
    :param minbinsize: The minimum bin size (default: 200,000).
    :return: None
    """
    outputdir = resultfile + '.filtersmallbins_' + str(minbinsize) + '.tsv'

    logger.info("Processing file:\t{}".format(fastafile))
    sequences = {}
    with open(fastafile, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if " " in line:
                    seq, others = line.split(' ', 1)
                    sequences[seq] = ""
                else:
                    seq = line.rstrip("\n")
                    sequences[seq] = ""
            else:
                sequences[seq] += line.rstrip("\n")
    logger.info("Reading Map:\t{}".format(resultfile))
    dic = {}
    bin_size = {}
    with open(resultfile, "r") as f:
        for line in f:
            contig_name, cluster_name = line.strip().split('\t')
            try:
                dic[cluster_name].append(contig_name)
                bin_size[cluster_name] += len(sequences['>' + contig_name])
            except:
                dic[cluster_name] = []
                dic[cluster_name].append(contig_name)
                bin_size[cluster_name] = len(sequences['>' + contig_name])

    with open(outputdir, "w") as f:
        for cluster_name in dic:
            if bin_size[cluster_name] >= minbinsize:
                for contigIdx in dic[cluster_name]:
                    f.write(contigIdx + "\t" + cluster_name + "\n")
    f.close()

    gen_bins(fastafile, outputdir, args.output_path + '/comebin_res_bins')
    shutil.copy2(outputdir, args.output_path + '/comebin_res.tsv')
