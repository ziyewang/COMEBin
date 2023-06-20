temperature=0.15
n_views=6
emb_szs_forcov=2048
emb_szs=2048
batch_size=1024
num_threads=48
output_path=??
out_augdata_path=??
bam_file_path=???
CUDA_VISIBLE_DEVICES=0

python main.py generate_aug_data --contig_file ${contig_file} \
--out_augdata_path ${out_augdata_path} \
--n_views ${n_views} --bam_file_path ${bam_file_path} --num_threads ${num_threads}


python main.py train --data ${out_augdata_path} \
--temperature ${temperature} --emb_szs_forcov ${emb_szs_forcov} \
--batch_size ${batch_size} --emb_szs ${emb_szs} --n_views ${n_views} \
--add_model_for_coverage \
--output_path ${output_path} --earlystop --addvars --vars_sqrt


#### (3) Clustering (run Leiden-based clustering methods and get the final result)

emb_file=${output_path}/embeddings.tsv
seed_file=${contig_file}.bacar_marker.2quarter_lencutoff_1001.seed

python main.py bin --contig_file ${contig_file} \
--emb_file ${emb_file} \
--output_path ${output_path} \
--seed_file ${seed_file} --num_threads ${num_threads}

python main.py get_result --contig_file ${contig_file} \
--output_path ${output_path} \
--seed_file ${seed_file} --num_threads ${num_threads}

