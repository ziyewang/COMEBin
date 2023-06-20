#####import
from Bio import SeqIO
import mimetypes
import os
import gzip
import random
import shutil

def get_inputsequences(fastx_file):
    file_type = mimetypes.guess_type(fastx_file)[1]
    if file_type == 'gzip':
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError("Unknown type of file: '{}".format(fastx_file))
    seqs = {}
    if os.path.getsize(fastx_file) == 0:
        return seqs
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
        seqs[seq_record.id] = seq_record.seq

    f.close()
    return seqs


def gen_augfasta(seqs, augprefix, out_file, p=None, contig_len=1000):
    seqkeys = []
    for seqid in seqs.keys():
        if len(seqs[seqid]) >= contig_len + 1:
            seqkeys.append(seqid)

    aug_seq_info = []
    if not p:
        with open(out_file, 'w') as f:
            for seqid in seqkeys:
                start = random.randint(0, len(seqs[seqid]) - (contig_len+1))
                sim_len = random.randint(contig_len, len(seqs[seqid]) - start)
                end = start + sim_len - 1
                # gen_seqs_dict[genome_name+"_sim_"+str(sim_count)] =seqs[seqid][start:end+1]
                sequence = str(seqs[seqid][start:end + 1])
                seqid_name = ">" + seqid + "_" + str(augprefix)
                f.write(seqid_name + "\n")
                f.write(sequence + "\n")
                aug_seq_info.append((seqid, start, end, sim_len))
    else:
        with open(out_file, 'w') as f:
            for seqid in seqkeys:
                sim_len = int(p * len(seqs[seqid]))
                start = random.randint(0, len(seqs[seqid]) - sim_len - 10)
                end = start + sim_len - 1
                # gen_seqs_dict[genome_name+"_sim_"+str(sim_count)] =seqs[seqid][start:end+1]
                sequence = str(seqs[seqid][start:end + 1])
                seqid_name = ">" + seqid + "_aug_" + str(augprefix)
                f.write(seqid_name + "\n")
                f.write(sequence + "\n")
                aug_seq_info.append((seqid, start, end, sim_len))

    aug_seq_info_out_file = out_file + '.aug_seq_info.tsv'

    with open(aug_seq_info_out_file, 'w') as afile:
        afile.write('seqid\tstart\tend\tlength\n')
        for i in range(len(aug_seq_info)):
            afile.write(
                aug_seq_info[i][0] + '\t' + str(aug_seq_info[i][1]) + '\t' + str(aug_seq_info[i][2]) + '\t' + str(
                    aug_seq_info[i][3]) + '\n')


########
def run_gen_augfasta(logger, args):
    # generate augmentation fasta file and save index
    num_aug = args.n_views - 1  # Generate several copies of augmented data
    fasta_file = args.contig_file
    out_path = args.out_augdata_path
    contig_len = args.contig_len

    outdir = out_path + '/aug0'
    os.makedirs(outdir)
    out_file = outdir + '/sequences_aug0.fasta'
    shutil.copyfile(fasta_file, out_file)

    from .gen_kmer import run_gen_kmer
    run_gen_kmer(out_file, 0, 4)

    for i in range(num_aug):
        outdir = out_path + '/aug' + str(i + 1)
        os.makedirs(outdir)
        logger.info("aug:\t" + str(i+1))
        p = None
        seqs = get_inputsequences(fasta_file)

        out_file = outdir + '/sequences_aug' + str(i + 1) + '.fasta'
        gen_augfasta(seqs, 'aug' + str(i + 1), out_file, p=p, contig_len=contig_len)
        run_gen_kmer(out_file, 0, 4)
