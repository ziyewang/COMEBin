# COMEBin
GitHub repository for the manuscript "COMEBin allows effective binning of metagenomic contigs using COntrastive Multi-viEw representation learning". 

## <a name="started"></a>Getting Started

### <a name="docker"></a>Install COMEBin via source code


Obtain codes and create an environment:
After installing Anaconda (or miniconda), first obtain COMEBin:

```sh
git clone https://github.com/ziyewang/COMEBin.git
```
Then create an environment to run COMEBin.

```sh
cd MetaBinner
conda env create -f comebin_env.yaml
conda activate comebin_env
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate bam files as input to our program.

Several binning methods can generate bam files by aligning reads to contigs (such as MetaWRAP and MetaBAT), and we provide one way to generate the input files as follows.
### Generate bam files 
To generate bam files from sequencing reads directly, run the following script slightly modified from the "binning.sh" of MetaWRAP. The script supports different types of sequencing reads, and the default type is "paired" ([readsX_1.fastq readsX_2.fastq ...]). 

```sh
cd path_to_COMEBin
cd COMEBin/scripts

bash gen_cov_file.sh -a contig_file \
-o output_dir_of_bamfiles \
path_to_sequencing_reads/*fastq

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory (to save the coverage files)
	-b STR          directory for the bam files (optional)
        -t INT          number of threads (default=1)
        -m INT          amount of RAM available (default=4)
        -l INT          minimum contig length to bin (default=1000bp).
        --single-end    non-paired reads mode (provide *.fastq files)
        --interleaved   input read files contain interleaved paired-end reads
        -f              Forward read suffix for paired reads (default="_1.fastq")
	-r              Reverse read suffix for paired reads (default="_2.fastq")

```

And the users can run the following command to keep the contigs longer than 1000bp for binning.

```sh
cd path_to_COMEBin
cd COMEBin/scripts

python Filter_tooshort.py final.contigs.fa 1000
```


## <a name="started"></a>An example to run COMEBin:

We ran COMEBin mainly in three steps: (a) Get augmentation data, (b) Get representation, and (c) Clustering (run Leiden-based clustering methods and get final result).
### (a) Get augmentation data,
```sh
time python main.py generate_aug_data --contig_file ${contig_file} \
--out_augdata_path ${out_augdata_path} \
--n_views 6 --bam_file_path ${bam_file_path} --num_threads 48
```
where ${bam_file_path} denotes the path to access the bam files



####################################################################
NokmerMetric
####################################################################
data=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/data_augmentation_clean
dataset_name=SRR5720233_4mer_6_view
output_path=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric

nepochs=200
temperature=0.15
n_views=6

#基本固定的参数
emb_szs_forcov=2048
emb_szs=2048
batch_size=1024
kmer=4mer
n_layer=3

time CUDA_VISIBLE_DEVICES=0 python main.py train --data ${data} \
--epochs ${nepochs} --temperature ${temperature} --emb_szs_forcov ${emb_szs_forcov} --dataset_name ${dataset_name} \
--kmer ${kmer} --batch_size ${batch_size} --emb_szs ${emb_szs} --n_views ${n_views} --n_layer ${n_layer} \
--add_model_for_coverage \
--output_path ${output_path} --earlystop --addvars --vars_sqrt

####################################################################
####binning
####################################################################
output_path=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric
emb_file=${output_path}/embeddings.tsv
contig_file=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/BATS_SAMN07137079_METAG.scaffolds.min500.fasta.f1k.fasta
seed_file=${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed

time python main.py bin --contig_file ${contig_file} \
--emb_file ${emb_file} \
--output_path ${output_path} \
--seed_file ${seed_file} --num_threads 48

#####
####################################################################
##gen_final_result
####################################################################
output_path=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric
emb_file=${output_path}/embeddings.tsv
contig_file=/mnt/data1/DeepBin/data/BATS_10samples/running_time/single_sample_mode_10sample/SRR5720233/BATS_SAMN07137079_METAG.scaffolds.min500.fasta.f1k.fasta
seed_file=${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed

time python main.py get_result --contig_file ${contig_file} \
--output_path ${output_path} \
--seed_file ${seed_file} --num_threads 48



```sh
```

## <a name="contact"></a>Contacts and bug reports
Please feel free to send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Prof. Shanfeng Zhu: zhusf@fudan.edu.cn


## <a name="References"></a>References
[3] https://github.com/dparks1134/UniteM.

[4] Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome research 25.7 (2015): 1043-1055.

[5] Christian M. K. Sieber, Alexander J. Probst., et al. (2018). "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy". Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1.

[6] Uritskiy, Gherman V., Jocelyne DiRuggiero, and James Taylor. "MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis." Microbiome 6.1 (2018): 1-13.

