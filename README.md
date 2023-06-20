# COMEBin
GitHub repository for the manuscript "MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities". We are glad that overall Metabinner achieve top performance on the CAMI II Challenge. Please refer to Meyer, F., et al[1] for the results of CAMI II Challenge.

MetaBinner consists of two modules: 1) “Component module” includes steps 1-4, developed for generating high-quality, diverse component binning results; and 2) “Ensemble module” includes step 5, developed for recovering individual genomes from the component binning results. MetaBinner is an ensemble binning method, but it does not need the outputs of other individual binners. Instead, MetaBinner generates multiple high-quality component binning results based on the proposed “partial seed” method, for further integration. Please see our manuscript for the details.

## <a name="started"></a>Getting Started

### <a name="docker"></a>Install COMEBin via source code


Obtain codes and create an environment:
After installing Anaconda (or miniconda), fisrt obtain COMEBin:

```sh
git clone https://github.com/ziyewang/COMEBin.git
```
Then simply create a environment to run COMEBin.

```sh
cd MetaBinner
conda env create -f comebin_env.yaml
conda activate comebin_env
```

## <a name="preprocessing"></a>Preprocessing

The preprocessing steps aim to generate bam files as input to our program.

There are several binning methods that can generate bam files by aligning reads to contigs (such as MetaWRAP and MetaBAT) and we provide one method to generate the input files as follows.
### Generate bam files 
To generate bam files from sequencing reads directly, run the following script slightly modified from the "binning.sh" of MetaWRAP. The script support different types of sequencing reads, and the defalut type is "paired" ([readsX_1.fastq readsX_2.fastq ...]). If MetaBinner is installed via bioconda, users can obtain path_to_MetaBinner via running this command: $(dirname $(which run_metabinner.sh))

```sh
cd path_to_COMEBin
cd scripts

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
        --interleaved   the input read files contain interleaved paired-end reads
        -f              Forward read suffix for paired reads (default="_1.fastq")
	-r              Reverse read suffix for paired reads (default="_2.fastq")

```

And the users can run the following command to keep the contigs longer than 1000bp for binning.

```
cd path_to_MetaBinner
cd scripts

python Filter_tooshort.py test_data/final.contigs_f1k.fa 1000
```


## <a name="started"></a>An example to run COMEBin:
```sh

#path to COMEBin
metabinner_path=/home/wzy/COMEBin
Note: If users intall MetaBinner via bioconda, they can set metabinner_path as follows: metabinner_path=$(dirname $(which run_metabinner.sh))

##test data
#path to the input files for MetaBinner and the output dir:
contig_file=/home/wzy/MetaBinner/test_data/final_contigs_f1k.fa
output_dir=/home/wzy/MetaBinner/test_data/output
coverage_profiles=/home/wzy/MetaBinner/test_data/coverage_profile_f1k.tsv
kmer_profile=/home/wzy/MetaBinner/test_data/kmer_4_f1000.csv


bash run_metabinner.sh -a ${contig_file} -o ${output_dir} -d ${coverage_profiles} -k ${kmer_profile} -p ${metabinner_path}

Options:

        -a STR          metagenomic assembly file
        -o STR          output directory
        -d STR          coverage_profile.tsv; The coverage profiles, containing a table where each row correspond
                            to a contig, and each column correspond to a sample. All values are separated with tabs.
        -k STR          kmer_profile.csv; The composition profiles, containing a table where each row correspond to a contig,
                            and each column correspond to the kmer composition of particular kmer. All values are separated with comma.
        -p STR          path to MetaBinner; e.g. /home/wzy/MetaBinner
        -t INT          number of threads (default=1)
        -s STR          Dataset scale; eg. small,large,huge (default:large); Users can choose "huge" to run MetaBinner on huge datasets
                        with lower memory requirements.

#The file "metabinner_result.tsv" in the "${output_dir}/metabinner_res" is the final output.
Note: all paths in the run_metabinner.sh options should be absolute
```

## <a name="contact"></a>Contacts and bug reports
Please feel free to send bug reports or questions to
Ziye Wang: zwang17@fudan.edu.cn and Prof. Shanfeng Zhu: zhusf@fudan.edu.cn


## <a name="References"></a>References
[1] Meyer, F., Fritz, A., Deng, ZL. et al. Critical Assessment of Metagenome Interpretation: the second round of challenges. Nat Methods (2022). https://doi.org/10.1038/s41592-022-01431-4

[3] https://github.com/dparks1134/UniteM.

[4] Parks, Donovan H., et al. "CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." Genome research 25.7 (2015): 1043-1055.

[5] Christian M. K. Sieber, Alexander J. Probst., et al. (2018). "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy". Nature Microbiology. https://doi.org/10.1038/s41564-018-0171-1.

[6] Uritskiy, Gherman V., Jocelyne DiRuggiero, and James Taylor. "MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis." Microbiome 6.1 (2018): 1-13.

