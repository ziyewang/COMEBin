#!/usr/bin/env python
import argparse
import gzip
import os


def gen_bins(fastafile, resultfile, outputdir):
    # read fasta file
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


# def gen_bins(fastafile, resultfile, outputdir, minbinsize=200000):
# read fasta file
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/gold_assembly/input/marmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/megahit_assembly/input/marmgCAMI2_short_read_pooled_megahit_assembly_f1k_clean.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_covembedding/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.015_nepoch200_earlystop_noedge75_addvars_vars_sqrt_covembeddings/cluster_res/test3_parallel_leiden_markerseed_RBConfigurationVertexPartition_nolength_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

#
# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='//home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group3/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group3/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_70_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011172/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011172/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011113/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011113/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011223/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011223/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011325/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011325/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011284/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011284/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011152/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011152/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011295/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011295/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011132/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011132/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011113/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011113/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011101/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011120/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#
#
# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_90_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_90_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/gold_assembly/input/marmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/gold_assembly/input/marmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/gold_assembly/input/marmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_50_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/megahit_assembly/input/marmgCAMI2_short_read_pooled_megahit_assembly_f1k_clean.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/megahit_assembly/input/marmgCAMI2_short_read_pooled_megahit_assembly_f1k_clean.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_110_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_70_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_20.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_addcovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_addcovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_addcovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_110_partgraph_ratio_20.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_110_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_90_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_90_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group3/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group3/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_50_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group3/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group3/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group3/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group3/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/vamb_m200000_rerun_0403/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/vamb_m200000_rerun_0403/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/vamb_m200000/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/vamb_m200000/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/vamb_m200000/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/vamb_m200000/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_nokmermetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group3/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group3/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/megahit_assembly/input/marmgCAMI2_short_read_pooled_megahit_assembly_f1k_clean.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/gold_assembly/input/marmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011132/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011132/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011152/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011152/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011172/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011172/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011284/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011284/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011223/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011223/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011295/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011295/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011325/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011325/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_covembedding/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_kmer_l2_normalize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_90_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_kmer_l2_normalize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_kmer_l2_normalize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_3layer2048/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_1layer1024/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_3layer2048/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_1layer1024/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_kmer_l2_normalize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_1layer1024/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_1layer1024/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.1_nepoch200_earlystop_addvars_noedge75_vars_sqrt_pretrain_not_load_kmermetric_state_finetunelr_ratio_1_kmermodel_1layer2048/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_8view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_kmer_l2_normalize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_90_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_4view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_8view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_8view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_70_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#
# # #
# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011113/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011113/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_50_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011101/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011101/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011120/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011120/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011132/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011132/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011152/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011152/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011172/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011172/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011223/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011223/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011284/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011284/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011295/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011295/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/single_sample/ERR011325/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/multi_sample_10sample/ERR011325/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/coassembly_10sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_4view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group2/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group2/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric_4view/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#
# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.07/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.12/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.15/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/megahit_assembly/input/rhimgCAMI2_short_read_pooled_megahit_f1k.clean.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.15/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt_Nokmermodel_addcovloss_covmodel_temperature0.07/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/plant_associated_dataset/gold_assembly/input/rhimgCAMI2_short_read_pooled_gsa_f1k.fasta'
# resultfile='/home/wzy/data/plant_associated_dataset/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.07/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/input/strmgCAMI2_short_read_pooled_megahit_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.15_covembeddings/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/input/strmgCAMI2_short_read_pooled_gold_standard_assembly_f1k.fasta'
# resultfile='/home/wzy/data/camistrainmadness_shortread/gold_assembly/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.07_covembedding/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt_Nokmermodel_addcovloss_covmodel_temperature0.07/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#

# minbinsize=200000
# fastafile='/home/wzy/data/camimarine_shortread/megahit_assembly/input/marmgCAMI2_short_read_pooled_megahit_assembly_f1k_clean.fasta'
# resultfile='/home/wzy/data/camimarine_shortread/megahit_assembly/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.15/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.12_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.12/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_noedge75_vars_sqrt_Nokmermodel_addcovloss_covmodel_temperature0.15_covembeddings/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_noedge75_vars_sqrt_Nokmermodel_addcovloss_covmodel_temperature0.07_covembeddings/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/MetaHIT_reads/input/final_assembly_f1k.clean.fasta'
# resultfile='/home/wzy/data/MetaHIT_reads/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_noedge75_addvars_vars_sqrt_NokmerMetric_addcovloss_covmodel_temperature0.15/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.3_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/vamb_m200000/COMEBIN_with_VAMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_5_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'
# resultfile='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_50_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_110_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_110_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_70_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_110_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_90_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_90_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_90_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'


# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nolength/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nolength_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker_nolength/res/test3_parallel_leiden_RBERVertexPartition_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_70_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/COMEBin_nocovloss_tau0.07_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res_nomarker/res/test3_parallel_leiden_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/clmb_m200000/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'
#
# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/clmb_m200000/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/clmb_m200000/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/res_from158/clmb_m200000/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/clmb_m200000/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/clmb_m200000_double_check/COMEBIN_with_CLMB_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges100respara_10_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
# resultfile='/home/wzy/data/cami_gt/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.2_res_maxedges100respara_1_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/STEC_data/for_new_method/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/STEC_data/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_5_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_skin/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_skin/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_10_partgraph_ratio_80.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_30_partgraph_ratio_100.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_mousegut/data_aug/data_augmentation_clean/aug0/sequences_aug0.fasta'
# resultfile='/home/wzy/data/cami_mousegut/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

# minbinsize=200000
# fastafile='/home/wzy/data/cami_10sample/input/gsa_f1k.fa'
# resultfile='/home/wzy/data/cami_10sample/output/semibin_Easy_coassemly_mode_v1.0.0_default/COMEBIN_with_semibin_embeddings_l2normaize/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_1_partgraph_ratio_50.tsv'
# outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'

minbinsize=200000
fastafile='/home/wzy/data/water_research_PRJNA542960/group1/data_augmentation_clean/aug0/sequences_aug0.fasta'
resultfile='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_nocovloss_tau0.15_nepoch200_earlystop_addvars_nedge75_vars_sqrt_NokmerMetric/cluster_res/res/test3_parallel_leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.05_res_maxedges100respara_50_partgraph_ratio_100.tsv'
outputdir=resultfile+'.filtersmallbins_'+str(minbinsize)+'.tsv'



print("Processing file:\t{}".format(fastafile))
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
print("Reading Map:\t{}".format(resultfile))
dic = {}
bin_size = {}
with open(resultfile, "r") as f:
    for line in f:
        contig_name, cluster_name = line.strip().split('\t')
        try:
            dic[cluster_name].append(contig_name)
            bin_size[cluster_name] += len(sequences['>'+contig_name])
        except:
            dic[cluster_name] = []
            dic[cluster_name].append(contig_name)
            bin_size[cluster_name] = len(sequences['>'+contig_name])

with open(outputdir, "w") as f:
    for cluster_name in dic:
        if bin_size[cluster_name] >= minbinsize:
            for contigIdx in dic[cluster_name]:
                f.write(contigIdx + "\t" + cluster_name + "\n")
f.close()


gen_bins(fastafile, outputdir, outputdir+'_bins')

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-f", help="original fasta file")
#     parser.add_argument("-r", help="tsv version result file")
#     parser.add_argument("-o", help="output dir")
#     args = parser.parse_args()
#     gen_bins(args.f, args.r, args.o)
