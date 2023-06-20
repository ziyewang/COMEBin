#!/usr/bin/python
# -*- coding:utf8 -*-
import sys,os,re
import click
from Bio import SeqIO

def process_file(reader):
    contigs = {}
    for seq_record in SeqIO.parse(reader, "fasta"):
        contigs[seq_record.id] = seq_record.seq
    return contigs

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('k', type=click.INT)
def main(input_file, k):
    with open(input_file) as fp:
        items = process_file(input_file)
    output_file = os.path.splitext(input_file)[0] + '_' + str(k) + '.fa'
    output_len_file = os.path.join(os.path.dirname(input_file), 'contig_length_filter' + str(k) + '.txt')
    with open(output_file, 'w') as fout1, open(output_len_file, 'w') as fout2:
        for key in items:
            length=int(len(items[key]))
            if length > k :
                fout1.write('>'+key+'\n'+str(items[key])+'\n')
                fout2.write(key+'\t'+str(length)+'\n')
    print('finished')

if __name__ == "__main__":
    #os.chdir("")
    main()
    
