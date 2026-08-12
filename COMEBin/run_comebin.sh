#!/usr/bin/env bash

##############################################################################################################################################################
# This script is meant to be run COMEBin after obtaining the bam files.
# Author of pipeline: Ziye Wang.
# For questions, bugs, and suggestions, contact me at zwang17@fudan.edu.cn
##############################################################################################################################################################
VERSION="1.1.0"

help_message () {
  echo ""
  echo "COMEBin version: $VERSION"
  echo "Usage: bash run_comebin.sh [options] -a contig_file -o output_dir -p bam_file_path"
	echo "Options:"
	echo ""
	echo "  -h              display this help message and exit"
	echo "  -V              display the COMEBin version and exit"
	echo "  -a STR          metagenomic assembly file"
	echo "  -o STR          output directory"
	echo "  -p STR          path to access to the bam files"
	echo "  -n INT          number of views for contrastive multiple-view learning (default=6)"
	echo "  -t INT          stage worker/thread count where explicitly supported (default=5)"
	echo "  -w INT          concurrent Leiden worker processes (default=same as -t)"
	echo "  -m INT          HNSW neighbors retained per contig (default=100)"
	echo "  -d STR          NN training device: cuda, cuda:N, or cpu (default=cuda)"
	echo "  -s INT          random seed for reproducible runs (default=not set)"
	echo "  -E FLOAT        HMMER sequence E-value used when marker HMMs do not define TC cutoffs (default=1e-5)"
	echo "  -l FLOAT        temperature in loss function (default=0.07 for assemblies with an N50 > 10000, default=0.15 for others)"
	echo "  -e INT          embedding size for comebin network (default=2048)"
	echo "  -c INT          embedding size for coverage network (default=2048)"
	echo "  -b INT          batch size for training process (default=1024)"
	echo "";}

report_stage_failure () {
    local stage_name=$1
    local stage_status=$2

    if (( stage_status >= 128 )); then
        echo "${stage_name} exited with status ${stage_status} (possible signal $((stage_status - 128)))." >&2
    else
        echo "${stage_name} exited with status ${stage_status}." >&2
    fi
    exit "${stage_status}"
}


########################################################################################################
########################     LOADING IN THE PARAMETERS AND RUNNING              ########################
########################################################################################################

num_threads=5
n_views=6
temperature=
emb_szs_forcov=2048
emb_szs=2048
batch_size=1024
device=cuda
random_seed=
leiden_workers=
max_edges=100
hmm_evalue=1e-5

while getopts hVa:o:p:n:t:w:m:d:s:E:l:e:c:b: OPT; do
 case ${OPT} in
  h) help_message
    exit 0
    ;;
  V) echo "COMEBin version: $VERSION"
    exit 0
    ;;
  a) contig_file=$(realpath -e -- "${OPTARG}") || exit 1
    ;;
  o) output_dir=$(realpath -m -- "${OPTARG}") || exit 1
    ;;
  p) bam_file_path=$(realpath -e -- "${OPTARG}") || exit 1
    ;;
  n) n_views=${OPTARG}
    ;;
  t) num_threads=${OPTARG}
    ;;
  w) leiden_workers=${OPTARG}
    ;;
  m) max_edges=${OPTARG}
    ;;
  d) device=${OPTARG}
    ;;
  s) random_seed=${OPTARG}
    ;;
  E) hmm_evalue=${OPTARG}
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

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
seed_options=()
if [ -n "${random_seed}" ]; then
    export PYTHONHASHSEED="${random_seed}"
    export CUBLAS_WORKSPACE_CONFIG="${CUBLAS_WORKSPACE_CONFIG:-:4096:8}"
    seed_options=(--seed "${random_seed}")
    echo "Random seed: ${random_seed}"
fi
echo "Leiden workers: ${leiden_workers:-same as -t (${num_threads})}"
echo "HNSW neighbors per contig: ${max_edges}"
echo "HMMER sequence E-value: ${hmm_evalue}"

# check parameter
if [ -z "${contig_file}" -o -z "${output_dir}" -o -z "${bam_file_path}" ]; then
  help_message
  exit 1
fi

script_path=$(realpath -- "${BASH_SOURCE[0]}") || exit 1
script_dir=$(dirname -- "${script_path}")

if [[ -f "${script_dir}/main.py" ]]; then
    comebin_dir="${script_dir}"
elif [[ -f "${script_dir}/COMEBin/main.py" ]]; then
    comebin_dir="${script_dir}/COMEBin"
else
    echo "Cannot locate COMEBin/main.py relative to ${script_path}." >&2
    exit 1
fi

main_py="${comebin_dir}/main.py"

temperature_options=()
if [ -n "${temperature}" ]; then
    temperature_options=(--temperature "${temperature}")
    echo "Tau(temperature): ${temperature}"
else
    echo "Tau(temperature): automatic from complete-assembly N50 in feature_manifest.json"
fi


########################################################################################################
###### Get augmentation data
########################################################################################################
folder="${output_dir}/data_augmentation"
augmentation_status=0
feature_manifest="${folder}/feature_manifest.json"
coverage_manifest="${folder}/coverage_manifest.json"

if [[ -s "${feature_manifest}" && -s "${folder}/contig_ids.txt" &&
      -s "${folder}/contig_lengths.npy" && -s "${folder}/augmentation_intervals.npy" &&
      -s "${folder}/kmer_counts.npy" && -s "${coverage_manifest}" &&
      -s "${folder}/bam_ids.txt" && -s "${folder}/coverage_mean.npy" &&
      -s "${folder}/coverage_var.npy" ]]; then
    echo "No need to run data augmentation."
else
    echo "Running data augmentation."
    python "${main_py}" generate_aug_data --contig_file "${contig_file}" \
    --out_augdata_path "${output_dir}/data_augmentation" \
    --n_views "${n_views}" --bam_file_path "${bam_file_path}" --num_threads "${num_threads}" "${seed_options[@]}"
    augmentation_status=$?
fi

if [[ ${augmentation_status} -ne 0 ]]; then
    report_stage_failure "Data augmentation" "${augmentation_status}"
fi

########################################################################################################
###### Get representation (training process)
########################################################################################################
train_status=0
embedding_dir="${output_dir}/comebin_res"
embedding_manifest="${embedding_dir}/embedding_manifest.json"
embedding_ids="${embedding_dir}/embedding_ids.txt"
emb_file="${embedding_dir}/embeddings.npy"
covemb_file="${embedding_dir}/covembeddings.npy"
embedding_lengths="${embedding_dir}/embedding_lengths.npy"

if [[ -s "${embedding_manifest}" && -s "${embedding_ids}" &&
      -s "${emb_file}" && -s "${covemb_file}" &&
      -s "${embedding_lengths}" ]]; then
    echo "No need to run getting representation."
else
    echo "Running getting representation."
    python "${main_py}" train --data "${output_dir}/data_augmentation" \
    --emb_szs_forcov "${emb_szs_forcov}" \
    --batch_size "${batch_size}" --emb_szs "${emb_szs}" --n_views "${n_views}" \
    --output_path "${output_dir}/comebin_res" --earlystop --addvars --vars_sqrt \
    --device "${device}" --num_threads "${num_threads}" "${temperature_options[@]}" "${seed_options[@]}"
    train_status=$?
fi


if [[ ${train_status} -ne 0 ]]; then
    report_stage_failure "Training" "${train_status}"
fi

########################################################################################################
#### Clustering (run Leiden-based clustering methods and get the final result)
########################################################################################################
seed_file="${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed"
leiden_options=(--max_edges "${max_edges}")
if [ -n "${leiden_workers}" ]; then
    leiden_options+=(--leiden_workers "${leiden_workers}")
fi
if [ -n "${random_seed}" ]; then
    leiden_options+=(--seed "${random_seed}")
fi
python "${main_py}" bin --contig_file "${contig_file}" \
--emb_file "${emb_file}" \
--output_path "${output_dir}/comebin_res" \
--seed_file "${seed_file}" --num_threads "${num_threads}" --hmm_evalue "${hmm_evalue}" "${leiden_options[@]}"
bin_status=$?

if [[ ${bin_status} -ne 0 ]]; then
    report_stage_failure "Clustering" "${bin_status}"
fi

python "${main_py}" get_result --contig_file "${contig_file}" \
--output_path "${output_dir}/comebin_res" \
--emb_file "${emb_file}" \
--seed_file "${seed_file}" --num_threads "${num_threads}" --max_edges "${max_edges}" \
--hmm_evalue "${hmm_evalue}" "${seed_options[@]}"
result_status=$?

if [[ ${result_status} -ne 0 ]]; then
    report_stage_failure "Final-result selection" "${result_status}"
fi
