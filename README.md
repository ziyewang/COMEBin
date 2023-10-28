# COMEBin
GitHub repository for the manuscript "COMEBin allows effective binning of metagenomic contigs using COntrastive Multi-viEw representation learning".
- [Overview](#overview)
- [System Requirements](#requirements)
- [Install COMEBin via bioconda](#install)
- [Install COMEBin via source code](#started)
- [A test dataset to demo COMEBin](#demo)
- [Preprocessing](#preprocessing)
- [How to run COMEBin](#runcomebin)
- [References](#References)
- [Contacts and bug reports](#contact)
  
## <a name="overview"></a>Overview
The framework of COMEBin is shown in the following figure, which is mainly divided into the following steps: 1) Data augmentation: construct five sets of augmented data by randomly extracting subsequences from the original contigs, resulting in six views for each original contig; 2) Construct feature vector: construct nucleotide frequency feature and coverage feature for each contig (including the original sequences and augmented sequences); 3) Contrastive learning: obtain low-dimensional embedding representations suitable for binning with heterogeneous information based on multi-view contrastive learning, and 4) Clustering: generate binning results based on community division algorithm Leiden.

Among them, the network structure used in contrastive learning includes two parts: 1) "Coverage network": process coverage features and 2) "Combine network": integrate the k-mer features and the "Coverage network" to obtain representations containing heterogeneous information.

<p align="center">
<img src="https://github.com/ziyewang/COMEBin/blob/master/overview.png" width="550"/>
</p>

## <a name="requirements"></a>System Requirements
### Hardware requirements
COMEBin requires only a standard computer with enough RAM to support the in-memory operations.

### OS Requirements
COMEBin v1.0.0 is supported and tested in Linux systems.

## <a name="install"></a>Install COMEBin via bioconda
To run COMEBin with CPU only:
```sh
conda create -n COMEBin_env
conda activate COMEBin_env
conda install -c conda-forge -c bioconda comebin
```
To run COMEBin with GPU (which provides faster performance when a GPU is available), you should also install PyTorch with GPU support:
```sh
conda create -n COMEBin_env
conda activate COMEBin_env
conda install -c conda-forge -c bioconda comebin
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia -c conda-forge
```

## <a name="started"></a>Install COMEBin via source code
You can also install COMEBin from the source code. 
After installing Anaconda (or miniconda), first, obtain COMEBin:

```sh
git clone https://github.com/ziyewang/COMEBin.git
```
Then, create an environment to run COMEBin.

```sh
cd path_to_COMEBin
conda env create -f comebin_env.yaml
conda activate comebin_env
```

Add dependencies: You need to run the command to make the file executable.
```sh
cd path_to_COMEBin/auxiliary
chmod +x test_getmarker_2quarter.pl
```

## <a name="demo"></a>A test dataset to demo COMEBin
We provide a small dataset to demo and test the software. Test data is available at https://drive.google.com/file/d/1xWpN2z8JTaAzWW4TcOl0Lr4Y_x--Fs5s/view?usp=sharing
```sh
comebin_test_data/BATS_SAMN07137077_METAG.scaffolds.min500.fasta.f1k.fasta
comebin_test_data/bamfiles/SRR5720343.bam
```
Run COMEBin on the test dataset:
```sh
cd path_to_COMEBin/COMEBin

CUDA_VISIBLE_DEVICES=0 bash run_comebin.sh -a path_to_comebin_test_data/BATS_SAMN07137077_METAG.scaffolds.min500.fasta.f1k.fasta \
-p path_to_comebin_test_data/bamfiles \
-o path_to_comebin_test_data/run_comebin_test \
-n 6 \
-t 40
```
Excepted output is given in: path_to_comebin_test_data/run_comebin_test

```sh
Final result (bins): path_to_comebin_test_data/run_comebin_test/comebin_res/comebin_res_bins
Final result in tsv format: path_to_comebin_test_data/run_comebin_test/comebin_res/comebin_res.tsv
```


## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate bam files as input to our program.

Several binning methods can generate bam files by aligning reads to contigs (such as MetaWRAP and MetaBAT), and we provide one way to generate the input files as follows.
### Generate bam files
To generate bam files from sequencing reads directly, run the script slightly modified from the "binning.sh" of MetaWRAP. The script supports different types of sequencing reads, and the default type is "paired" ([readsX_1.fastq readsX_2.fastq ...]).

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
        -l INT          minimum contig length to the bin (default=1000bp).
        --single-end    non-paired reads mode (provide *.fastq files)
        --interleaved   input read files contain interleaved paired-end reads
        -f              Forward read suffix `for` paired reads (default="_1.fastq")
	-r              Reverse read suffix for paired reads (default="_2.fastq")

```

And the users can run the following command to keep the contigs longer than 1000bp for binning.

```sh
cd path_to_COMEBin
cd COMEBin/scripts

python Filter_tooshort.py final.contigs.fa 1000
```


## <a name="runcomebin"></a>How to run COMEBin
### Run COMEBin via bioconda
```sh
run_comebin.sh -a ${contig_file} \
-o ${output_path} \
-p ${path_to_bamfiles} \
-t 40
```

### Run COMEBin via source code
```sh
cd path_to_COMEBin/COMEBin

bash run_comebin.sh -a ${contig_file} \
-o ${output_path} \
-p ${path_to_bamfiles} \
-t 40
```
```sh
Usage: run_comebin.sh [options] -a contig_file -o output_dir -p bam_file_path
Options:

  -a STR          metagenomic assembly file
  -o STR          output directory
  -p STR          path to access to the bam files
  -n INT          number of views for contrastive multiple-view learning (default=6)
  -t INT          number of threads (default=5)
  -l FLOAT        temperature in loss function (default=0.07 for assemblies with an N50 > 10000, default=0.15 for others)
  -e INT          embedding size for combine network (default=2048)
  -c INT          embedding size for coverage network (default=2048)
  -b INT          batch size for training process (default=1024)
```

## <a name="References"></a>References
[1] Meyer, F., Fritz, A., Deng, ZL. et al. Critical Assessment of Metagenome Interpretation: the second round of challenges. Nature Methods (2022). https://doi.org/10.1038/s41592-022-01431-4

[2] Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome Research 25.7 (2015): 1043-1055.

[3] https://github.com/dparks1134/UniteM.

[4] Pan, S., Zhu, C., Zhao, XM., and Coelho, LP. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nature Communications 13, 2326 (2022). https://doi.org/10.1038/s41467-022-29843-y

## <a name="contact"></a>Contacts and bug reports
Please feel free to send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Prof. Shanfeng Zhu: zhusf@fudan.edu.cn


