#!/usr/bin/env bash

##############################################################################################################################################################
# This script is meant to be run COMEBin after obtaining the bam files.
# Author of pipeline: Ziye Wang.
# For questions, bugs, and suggestions, contact me at zwang17@fudan.edu.cn
##############################################################################################################################################################
VERSION="1.0.2"

help_message () {
  echo ""
  echo "COMEBin version: $VERSION"
  echo "Usage: bash run_comebin.sh [options] -a contig_file -o output_dir -p bam_file_path"
	echo "Options:"
	echo ""
	echo "  -a STR          metagenomic assembly file"
	echo "  -o STR          output directory"
	echo "  -p STR          path to access to the bam files"
	echo "  -n INT          number of views for contrastive multiple-view learning (default=6)"
	echo "  -t INT          number of threads (default=5)"
	echo "  -l FLOAT        temperature in loss function (default=0.07 for assemblies with an N50 > 10000, default=0.15 for others)"
	echo "  -e INT          embedding size for comebin network (default=2048)"
	echo "  -c INT          embedding size for coverage network (default=2048)"
	echo "  -b INT          batch size for training process (default=1024)"
	echo "";}


########################################################################################################
########################     LOADING IN THE PARAMETERS AND RUNNING              ########################
########################################################################################################

num_threads=5
n_views=6
#temperature=0.15
emb_szs_forcov=2048
emb_szs=2048
batch_size=1024

while getopts a:o:p:n:t:l:e:c:b: OPT; do
 case ${OPT} in
  a) contig_file=${OPTARG}
    ;;
  o) output_dir=${OPTARG}
    ;;
  p) bam_file_path=${OPTARG}
    ;;
  n) n_views=${OPTARG}
    ;;
  t) num_threads=${OPTARG}
    ;;
  l) temperature=${OPTARG}
    ;;
  e) emb_szs=${OPTARG}
    ;;
  c) emb_szs_forcov=${OPTARG}
    ;;
  b) batch_size=${OPTARG}
    ;;
  \?)
#    printf "[Usage] `date '+%F %T'` -i <INPUT_FILE> -o <OUTPUT_DIR> -o <P
#RODUCT_CODE> -s <SOFTWARE_VERSION> -t <TYPE>\n" >&2
    exit 1
 esac
done

# check parameter
if [ -z "${contig_file}" -o -z "${output_dir}" -o -z "${bam_file_path}" ]; then
  help_message
  exit 1
fi


if [ -z "$temperature" ]; then
    # Compute the length of each sequence and sort using the awk command
    awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$contig_file" | sort -rn > ${contig_file}_lengths.txt

    # CAL N50
    total_length=$(awk '{sum+=$1} END {print sum}' ${contig_file}_lengths.txt)
    target_length=$(awk -v total="$total_length" 'BEGIN {cutoff=total/2; current=0} {current+=$1; if (current >= cutoff) {print $1; exit}}' ${contig_file}_lengths.txt)

    # N50
    echo "N50: $target_length"
    # Check if N50 is greater than 10000 and set tau accordingly
    if [ "$target_length" -gt 10000 ]; then
        temperature=0.07
    else
        temperature=0.15
    fi
    echo "Tau(temperature): ${temperature}"
else
    echo "Tau(temperature): ${temperature}"
fi


current_path=$(pwd)
chmod +x ${current_path}/../auxiliary/test_getmarker_2quarter.pl

########################################################################################################
###### Get augmentation data
########################################################################################################
folder=${output_dir}/data_augmentation
keyword="_datacoverage_mean"

if [ -d "$folder" ]; then
    echo "${output_dir}/data_augmentation exists."
    count=$(find "$folder" -maxdepth 1 -type f -name "*$keyword*" | wc -l)
    echo "Number of files containing '$keyword' in the folder: $count"
    if [ "$count" -ne ${n_views} ]; then
        echo "Running data augmentation."
        python main.py generate_aug_data --contig_file ${contig_file} \
        --out_augdata_path ${output_dir}/data_augmentation \
        --n_views ${n_views} --bam_file_path ${bam_file_path} --num_threads ${num_threads}
    else
        echo "No need to run data augmentation."
    fi
else
    echo "${output_dir}/data_augmentation does not exist."
    echo "Running data augmentation."
    python main.py generate_aug_data --contig_file ${contig_file} \
    --out_augdata_path ${output_dir}/data_augmentation \
    --n_views ${n_views} --bam_file_path ${bam_file_path} --num_threads ${num_threads}
fi

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running generating augmentation data. Exiting.";exit 1; fi

########################################################################################################
###### Get representation (training process)
########################################################################################################
folder=${output_dir}/comebin_res
keyword="embeddings.tsv"

if [ -d "$folder" ]; then
    echo "${output_dir}/comebin_res exists."
    count=$(find "$folder" -maxdepth 1 -type f -name "*$keyword*" | wc -l)
    echo "Number of files containing '$keyword' in the folder: $count"
    if [ "$count" -ne 2 ]; then
        echo "Running getting representation."
        python main.py train --data ${output_dir}/data_augmentation \
        --temperature ${temperature} --emb_szs_forcov ${emb_szs_forcov} \
        --batch_size ${batch_size} --emb_szs ${emb_szs} --n_views ${n_views} \
        --add_model_for_coverage \
        --output_path ${output_dir}/comebin_res --earlystop --addvars --vars_sqrt
    else
        echo "No need to run getting representation."
    fi
else
    echo "${output_dir}/comebin_res does not exist."
    echo "Running getting representation."
    python main.py train --data ${output_dir}/data_augmentation \
    --temperature ${temperature} --emb_szs_forcov ${emb_szs_forcov} \
    --batch_size ${batch_size} --emb_szs ${emb_szs} --n_views ${n_views} \
    --add_model_for_coverage \
    --output_path ${output_dir}/comebin_res --earlystop --addvars --vars_sqrt
fi


if [[ $? -ne 0 ]] ; then echo "Something went wrong with running training network. Exiting.";exit 1; fi

########################################################################################################
#### Clustering (run Leiden-based clustering methods and get the final result)
########################################################################################################
emb_file=${output_dir}/comebin_res/embeddings.tsv
seed_file=${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed

python main.py bin --contig_file ${contig_file} \
--emb_file ${emb_file} \
--output_path ${output_dir}/comebin_res \
--seed_file ${seed_file} --num_threads ${num_threads}

python main.py get_result --contig_file ${contig_file} \
--output_path ${output_dir}/comebin_res \
--seed_file ${seed_file} --num_threads ${num_threads}

if [[ $? -ne 0 ]] ; then echo "Something went wrong with running clustering. Exiting.";exit 1; fi

