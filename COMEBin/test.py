sc.pp.neighbors(adata, n_neighbors=50, n_pcs=0, use_rep='X')
sc.tl.leiden(adata,resolution=50)
pred = adata.obs['leiden'].to_list()
pred = [int(x) for x in pred]


#
#STEC
output_path='/home/wzy/data/STEC_data/output/COMEBin_whole_addcovloss_batchsize_512/'
seed_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr16_17-07-50_ZzStudio-7048-4x1080/input/final_assembly_f1k.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
emb_file=output_path+'/embeddings.tsv'
contig_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr16_17-07-50_ZzStudio-7048-4x1080/input/final_assembly_f1k.fasta'

# ### STEC co-assembly
# output_path='/home/wzy/data/STEC_data/coassembly_5sample/output/COMEBin_whole_addcovloss_batchsize_512'
# seed_file='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
# emb_file=output_path+'/embeddings.tsv'
# contig_file='/home/wzy/data/STEC_data/coassembly_5sample/input/final_assembly.f1k.fasta'


# ####cami_gt
# output_path='/home/wzy/data/cami_gt/output/COMEBin_whole_addcovloss_batchsize_512/'
# seed_file='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
# emb_file=output_path+'embeddings.tsv'
# contig_file='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'
#
# # ##
# ####watergroup1
# output_path='/home/wzy/data/water_research_PRJNA542960/group1/output/COMEBin_whole_addcovloss_batchsize_512'
# seed_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/Apr23_02-12-40_ZzStudio-7048-4x1080/input/final.contigs.fa.f1k.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
# emb_file=output_path+'/embeddings.tsv'
# contig_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/Apr23_02-12-40_ZzStudio-7048-4x1080/input/final.contigs.fa.f1k.fasta'


# ##cami_mousegut
# output_path='/home/wzy/data/cami_mousegut/output/COMEBin_whole_addcovloss_batchsize_512/'
# seed_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr19_09-30-42_ZzStudio-7048-4x1080/input/anonymous_gsa_pooled_f1k.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
# emb_file=output_path+'embeddings.tsv'
# contig_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr19_09-30-42_ZzStudio-7048-4x1080/input/anonymous_gsa_pooled_f1k.fasta'


emb_file2=output_path+'covembeddings.tsv'


output_path = output_path + '/cluster_res/'

embHeader = pd.read_csv(emb_file, sep='\t', nrows=1)
embMat = pd.read_csv(emb_file, sep='\t', usecols=range(1, embHeader.shape[1])).values
namelist = pd.read_csv(emb_file, sep='\t', usecols=range(1)).values[:, 0]

lengths = get_length(contig_file)
length_weight = []
for seq_id in namelist:
    length_weight.append(lengths[seq_id])

norm_embeddings = normalize(embMat)

embHeader2 = pd.read_csv(emb_file2, sep='\t', nrows=1)
embMat2 = pd.read_csv(emb_file2, sep='\t', usecols=range(1, embHeader2.shape[1])).values
norm_embeddings2 = normalize(embMat2)
# norm_embeddings=norm_embeddings[:10]

# max_edges_list = list(np.arange(0,320,20))
# max_edges_list[0] = 10



def run_leiden(logger,norm_embeddings,max_edges_list,seed_file, namelist, length_weight=None, partition_type='SurpriseVertexPartition',resolution_parameter=None):
    # max_edges=100
    # max_edges_list=[10,30,50,70,90,110,130,150,170,200,220,250]

    for max_edges in max_edges_list:
        p = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10)


        time_start = time.time()
        ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges+1)
        #ann_distances是l2 distance的平方
        time_end = time.time()
        logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
        sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
        targets_indices = ann_neighbor_indices[:,1:]
        targets = targets_indices.flatten()
        wei = ann_distances[:,1:]

        wei_save = wei.copy()

        sigmas_sq = np.median(wei, axis=1)
        sigmas = np.sqrt(sigmas_sq)

        for i in range(len(wei)):
            num = 2 * sigmas[i] * sigmas[targets_indices[i]]
            den = sigmas_sq[i] + sigmas_sq[targets_indices[i]]
            wei_save[i] = np.sqrt(num/den) * np.exp(-wei[i]/ den)
            # print(wei[i])

        seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
        initial_list = list(np.arange(len(namelist)))
        is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

        wei_save = wei_save.flatten()

        index = sources > targets
        sources = sources[index]
        targets = targets[index]
        wei_save = wei_save[index]
        vcount = len(norm_embeddings)
        edgelist = list(zip(sources, targets))
        g = Graph(vcount, edgelist)

        if partition_type == 'SurpriseVertexPartition':
            res = leidenalg.SurpriseVertexPartition(g,
                                                    weights=wei_save, initial_membership = initial_list,
                                                    node_sizes=length_weight)
            nodesize = True
        if partition_type == 'RBConfigurationVertexPartition':
            res = leidenalg.RBConfigurationVertexPartition(g,
                                                           weights=wei_save, initial_membership = initial_list,
                                                           resolution_parameter = resolution_parameter)
            nodesize = False

        if partition_type == 'RBERVertexPartition':
            res = leidenalg.RBERVertexPartition(g,
                                                weights=wei_save, initial_membership = initial_list,
                                                resolution_parameter = resolution_parameter,node_sizes=length_weight)

            nodesize = True

        optimiser = leidenalg.Optimiser()
        diff = 1
        while diff > 0:
            diff = optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed,n_iterations=-1)


        part = list(res)

        modularity = res.modularity
        if partition_type=='SurpriseVertexPartition':
            output_file = output_path + 'leiden_markerseed_partition_type_'+partition_type+'_nodesize_'+str(nodesize)+'_max_edges_'+str(max_edges)+'_modularity_'+str(modularity)+'.tsv'
        else:
            output_file = output_path + 'leiden_markerseed_partition_type_'+partition_type+'_nodesize_'+str(nodesize)\
                          +'_max_edges_'+str(max_edges)+'_respara_'+str(resolution_parameter)+'_modularity_'+str(modularity)+'.tsv'

        contig_labels_dict ={}
        # dict of communities
        numnode = 0
        rang = []
        for ci in range(len(part)):
            rang.append(ci)
            numnode = numnode+len(part[ci])
            for id in part[ci]:
                contig_labels_dict[namelist[id]] = 'group'+str(ci)

        print(output_file)
        f = open(output_file, 'w')
        for contigIdx in range(len(contig_labels_dict)):
            f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
        f.close()


#####
## run RBConfigurationVertexPartition
max_edges_list = list(np.arange(0,300,50))
max_edges_list[0] = 10
parameter_list = list(np.arange(0,150,20))
parameter_list[0] = 5

for res_para in parameter_list:
    run_leiden(logger,norm_embeddings,max_edges_list,seed_file, namelist, length_weight=None,
               partition_type='RBConfigurationVertexPartition', resolution_parameter=res_para)

### run RBERVertexPartition nodesize
max_edges_list = list(np.arange(0,300,50))
max_edges_list[0] = 10
parameter_list = list(np.arange(0,150,20))
parameter_list[0] = 5

for res_para in parameter_list:
    run_leiden(logger,norm_embeddings,max_edges_list,seed_file, namelist, length_weight=length_weight,
               partition_type='RBERVertexPartition', resolution_parameter=res_para)


#https://github.com/scverse/scanpy/blob/d7e13025b931ad4afd03b4344ef5ff4a46f78b2b/scanpy/neighbors/__init__.py#L434

"""

from sklearn.neighbors import kneighbors_graph
from igraph import Graph
import os
import pandas as pd
import numpy as np
import gzip
from Bio import SeqIO
import mimetypes
from sklearn.preprocessing import normalize

import functools
from sklearn.cluster.k_means_ import euclidean_distances, stable_cumsum, KMeans, check_random_state, row_norms, MiniBatchKMeans
import scipy.sparse as sp

import hnswlib
import time
# import leidenalg

def get_graphmatrix(data, num_process=40):
    time_start=time.time()
    distance_matrix = kneighbors_graph(
        data,
        n_neighbors=min(max_edges, data.shape[0] - 1),
        mode='distance',
        p=2,
        n_jobs=num_process)
    time_end = time.time()
    print('time cost', time_end-time_start,"s")
    return distance_matrix

max_edges = 3
a= get_graphmatrix(norm_embeddings[:10], num_process=40)
Dsq=a.power(2)
W = (
    Dsq.copy()
)  # need to copy the distance matrix here; what follows is inplace
for i in range(len(Dsq.indptr[:-1])):
    row = Dsq.indices[Dsq.indptr[i] : Dsq.indptr[i + 1]]
    num = 2 * sigmas[i] * sigmas[row]
    den = sigmas_sq[i] + sigmas_sq[row]
    W.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] = np.sqrt(num / den) * np.exp(
        -Dsq.data[Dsq.indptr[i] : Dsq.indptr[i + 1]] / den
    )
    
"""
import scanpy as sc

# max_edges_list = list(np.arange(0,300,50))
# max_edges_list[0] = 10
# max_edges_list.append(5)
# max_edges_list.append(15)

max_edges_list = [100,150,200,350,300,350,400,50,]
# parameter_list = list(np.arange(0,150,20))
# parameter_list[0] = 5RBERVertexPartition
parameter_list=[160,180,200,220,240,260,280,300]

for max_edges in max_edges_list:
    for res_para in parameter_list:
        adata = sc.AnnData(norm_embeddings)
        sc.pp.neighbors(adata, n_neighbors=max_edges, n_pcs=0, use_rep='X')
        sc.tl.leiden(adata,resolution=res_para)
        pred = adata.obs['leiden'].to_list()
        pred = [int(x) for x in pred]
        output_file=output_path+'/test_leiden.scanpy_edge_'+str(max_edges)+'_resolu_'+str(res_para) +'_default.tsv'
        f = open(output_file, 'w')
        for contigIdx in range(len(pred)):
            f.write(namelist[contigIdx] + "\t" + str(pred[contigIdx]) + "\n")
        f.close()



######
#mu
max_edges = 200
p = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10, space='cosine')

def cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_multiplex(logger, length_weight, seed_file, max_edges, norm_embeddings, p, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    time_start = time.time()
    ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")

    partgraph_ratio_list =[99,80,50]



    graphs = []
    weis = []
    for partgraph_ratio in partgraph_ratio_list:
        sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
        targets_indices = ann_neighbor_indices[:,1:]
        targets = targets_indices.flatten()

        wei = ann_distances[:,1:]
        wei = wei.flatten()

        wei = 1 - wei
        dist_cutoff = np.percentile(wei, partgraph_ratio)
        print(dist_cutoff)
        save_index = wei >= dist_cutoff

        sources = sources[save_index]
        targets = targets[save_index]
        wei = wei[save_index]

        index = sources > targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        vcount = len(norm_embeddings)
        edgelist = list(zip(sources, targets))
        g = Graph(vcount, edgelist)
        graphs.append(g)
        weis.append(wei)

    res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                            weights=weis[0], initial_membership = initial_list,
                                            node_sizes=length_weight)
    res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                             weights=weis[1], initial_membership = initial_list,
                                             node_sizes=length_weight)
    res3 = leidenalg.SurpriseVertexPartition(graphs[2],
                                             weights=weis[2], initial_membership = initial_list,
                                             node_sizes=length_weight)


    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition_multiplex([res1, res2, res3], is_membership_fixed=is_membership_fixed,n_iterations=-1)


    part = list(res1)

    return part, bandwidth


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res1, is_membership_fixed=is_membership_fixed,n_iterations=-1)

part = list(res1)

output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_maxedge200_partgraph_1.tsv'

contig_labels_dict ={}
# dict of communities
numnode = 0
rang = []
for ci in range(len(part)):
    rang.append(ci)
    numnode = numnode+len(part[ci])
    for id in part[ci]:
        contig_labels_dict[namelist[id]] = 'group'+str(ci)

print(output_file)
f = open(output_file, 'w')
for contigIdx in range(len(contig_labels_dict)):
    f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
f.close()

###试试串联
res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                         weights=weis[0], initial_membership = initial_list,
                                         node_sizes=length_weight)

optimiser = leidenalg.Optimiser()
diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res1, is_membership_fixed=is_membership_fixed,n_iterations=-1)


res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                         weights=weis[1], initial_membership = res1.membership,
                                         node_sizes=length_weight)

optimiser = leidenalg.Optimiser()
diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res2, is_membership_fixed=is_membership_fixed,n_iterations=-1)

res3 = leidenalg.SurpriseVertexPartition(graphs[2],
                                         weights=weis[2], initial_membership = res2.membership,
                                         node_sizes=length_weight)

optimiser = leidenalg.Optimiser()
diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res3, is_membership_fixed=is_membership_fixed,n_iterations=-1)



part = list(res3)


#####

def cluster_hnswlib_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_multiplex(logger, length_weight, seed_file, max_edges, norm_embeddings, p, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l2', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=100 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    time_start = time.time()
    ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")

    bandwidth_list = [0.1,0.2,0.3]

    weis = []
    graphs = []

    for bandwidth in bandwidth_list:
        sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
        targets_indices = ann_neighbor_indices[:,1:]
        targets = targets_indices.flatten()
        wei = ann_distances[:,1:]
        wei = wei.flatten()
    
        dist_cutoff = np.percentile(wei, partgraph_ratio)
        save_index = wei <= dist_cutoff
        sources = sources[save_index]
        targets = targets[save_index]
        wei = wei[save_index]

        if lmode == 'l1':
            wei = np.sqrt(wei)
            if not bandwidth:
                bandwidth = np.percentile(wei, percentile)
            wei = np.exp(-wei / bandwidth) #updated

        if lmode == 'l2':
            if not bandwidth:
                bandwidth = np.percentile(wei, percentile)
            wei = np.exp(-wei / bandwidth)

        logger.info("bandwidth:\t" +str(bandwidth))

        index = sources > targets
        sources = sources[index]
        targets = targets[index]
        wei = wei[index]
        vcount = len(norm_embeddings)
        edgelist = list(zip(sources, targets))
        g = Graph(vcount, edgelist)

        weis.append(wei)
        graphs.append(g)

    res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                             weights=weis[0], initial_membership = initial_list,
                                             node_sizes=length_weight)
    res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                             weights=weis[1], initial_membership = initial_list,
                                             node_sizes=length_weight)
    res3 = leidenalg.SurpriseVertexPartition(graphs[2],
                                             weights=weis[2], initial_membership = initial_list,
                                             node_sizes=length_weight)


    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition_multiplex([res1, res2, res3], is_membership_fixed=is_membership_fixed,n_iterations=-1)

    part = list(res1)

    return part, bandwidth

cluster_hnswlib_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_multiplex(logger, length_weight, seed_file,
                                                                                       max_edges, norm_embeddings,
                                                                                       p, namelist, num_process=40, percentile=5,
                                                                                       bandwidth=None, lmode='l1', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=100)


res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                         weights=weis[0], initial_membership = initial_list,
                                         node_sizes=length_weight)


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition_multiplex(res1, is_membership_fixed=is_membership_fixed,n_iterations=-1)

part = list(res1)

res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                         weights=weis[0], initial_membership = initial_list,
                                         node_sizes=length_weight)


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res1, is_membership_fixed=is_membership_fixed,n_iterations=-1)

part = list(res1)

res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                         weights=weis[1], initial_membership = initial_list,
                                         node_sizes=length_weight)


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res2, is_membership_fixed=is_membership_fixed,n_iterations=-1)

part = list(res2)

# output_file = output_path+'cluster_hnswlib_l1mode_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_multiplex.tsv'

# output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_emb_covemb_maxedege200_partgraph_100.tsv'
output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_covemb_maxedege200_partgraph_100.tsv'


contig_labels_dict ={}
# dict of communities
numnode = 0
rang = []
for ci in range(len(part)):
    rang.append(ci)
    numnode = numnode+len(part[ci])
    for id in part[ci]:
        contig_labels_dict[namelist[id]] = 'group'+str(ci)

print(output_file)
f = open(output_file, 'w')
for contigIdx in range(len(contig_labels_dict)):
    f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
f.close()



####
max_edges = 200
p1 = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10, space='cosine')
p2 = fit_hnsw_index(logger, norm_embeddings2, ef=max_edges * 10, space='cosine')

def cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_emb1_covemb_multiplex(logger, length_weight, seed_file, max_edges, norm_embeddings, norm_embeddings2,p1,p2, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    graphs = []
    weis = []

    time_start = time.time()
    ann_neighbor_indices1, ann_distances1 = p1.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources1 = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices1 = ann_neighbor_indices1[:,1:]
    targets1 = targets_indices1.flatten()


    wei1 = ann_distances1[:,1:]
    wei1 = wei1.flatten()

    wei1 = 1 - wei1
    dist_cutoff = np.percentile(wei1, partgraph_ratio)
    print(dist_cutoff)
    save_index = wei1 >= dist_cutoff

    sources1 = sources1[save_index]
    targets1 = targets1[save_index]
    wei1 = wei1[save_index]

    index1 = sources1 > targets1
    sources1 = sources1[index1]
    targets1 = targets1[index1]
    wei1 = wei1[index1]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources1, targets1))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei1)

    ####cov graph
    time_start = time.time()
    ann_neighbor_indices2, ann_distances2 = p2.knn_query(norm_embeddings2, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources2 = np.repeat(np.arange(len(norm_embeddings2)), max_edges)
    targets_indices2 = ann_neighbor_indices2[:,1:]
    targets2 = targets_indices2.flatten()

    wei2 = ann_distances2[:,1:]
    wei2 = wei2.flatten()

    wei2 = 1 - wei2
    dist_cutoff = np.percentile(wei2, partgraph_ratio)
    print(dist_cutoff)
    save_index = wei2 >= dist_cutoff

    sources2 = sources2[save_index]
    targets2 = targets2[save_index]
    wei2 = wei2[save_index]

    index2 = sources2 > targets2
    sources2 = sources2[index2]
    targets2 = targets2[index2]
    wei2 = wei2[index2]
    vcount = len(norm_embeddings2)
    edgelist = list(zip(sources2, targets2))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei2)

    res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                             weights=weis[0], initial_membership = initial_list,
                                             node_sizes=length_weight)
    res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                             weights=weis[1], initial_membership = initial_list,
                                             node_sizes=length_weight)

    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition_multiplex([res1, res2], is_membership_fixed=is_membership_fixed,
                                                     n_iterations=-1)


    part = list(res1)

    return part, bandwidth


def cluster_hnswlib_cosine_leiden_markerseed_RBERVertexPartition_partgraph_emb1_covemb_multiplex(logger, length_weight, seed_file, max_edges, norm_embeddings, norm_embeddings2,p1,p2, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    graphs = []
    weis = []

    time_start = time.time()
    ann_neighbor_indices1, ann_distances1 = p1.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources1 = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices1 = ann_neighbor_indices1[:,1:]
    targets1 = targets_indices1.flatten()


    wei1 = ann_distances1[:,1:]
    wei1 = wei1.flatten()

    wei1 = 1 - wei1
    dist_cutoff = np.percentile(wei1, partgraph_ratio)
    print(dist_cutoff)
    save_index = wei1 >= dist_cutoff

    sources1 = sources1[save_index]
    targets1 = targets1[save_index]
    wei1 = wei1[save_index]

    index1 = sources1 > targets1
    sources1 = sources1[index1]
    targets1 = targets1[index1]
    wei1 = wei1[index1]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources1, targets1))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei1)

    ####cov graph
    time_start = time.time()
    ann_neighbor_indices2, ann_distances2 = p2.knn_query(norm_embeddings2, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources2 = np.repeat(np.arange(len(norm_embeddings2)), max_edges)
    targets_indices2 = ann_neighbor_indices2[:,1:]
    targets2 = targets_indices2.flatten()

    wei2 = ann_distances2[:,1:]
    wei2 = wei2.flatten()

    wei2 = 1 - wei2
    dist_cutoff = np.percentile(wei2, partgraph_ratio)
    print(dist_cutoff)
    save_index = wei2 >= dist_cutoff

    sources2 = sources2[save_index]
    targets2 = targets2[save_index]
    wei2 = wei2[save_index]

    index2 = sources2 > targets2
    sources2 = sources2[index2]
    targets2 = targets2[index2]
    wei2 = wei2[index2]
    vcount = len(norm_embeddings2)
    edgelist = list(zip(sources2, targets2))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei2)

    res1 = leidenalg.RBERVertexPartition(graphs[0],
                                             weights=weis[0], initial_membership = initial_list,
                                         node_sizes=length_weight)
    res2 = leidenalg.RBERVertexPartition(graphs[1],
                                             weights=weis[1], initial_membership = initial_list,
                                         node_sizes=length_weight)

    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition_multiplex([res1, res2], is_membership_fixed=is_membership_fixed,
                                                      n_iterations=-1)


    part = list(res1)

    return part, bandwidth


output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_RBERVertexPartition_nodesize_partgraph_emb_covemb_maxedege200_partgraph_100.tsv'
output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_RBERVertexPartition_nodesize_partgraph_emb_maxedege200_partgraph_100.tsv'
output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_RBERVertexPartition_nodesize_partgraph_covemb_maxedege200_partgraph_100.tsv'

output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_covemb_emb_chuanlian_maxedege200_partgraph_50.tsv'
#output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_maxedege200_partgraph_50.tsv'
#output_file = output_path+'cluster_hnswlib_cosine_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_covemb_maxedege200_partgraph_50.tsv'

output_file = output_path+'leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.15_res_maxedges200respara_30_partgraph_ratio_20.tsv.chuanlian_res1nofixnodes.tsv'

output_file = output_path+'leiden_markerseed_RBERVertexPartition_nodesize_hnsw_l2mode_bandwidth_0.1_res_maxedges50respara_30.tsv.addHeader.tsv.chuanlian.tsv'

contig_labels_dict ={}
# dict of communities
numnode = 0
rang = []
for ci in range(len(part)):
    rang.append(ci)
    numnode = numnode+len(part[ci])
    for id in part[ci]:
        contig_labels_dict[namelist[id]] = 'group'+str(ci)

print(output_file)
f = open(output_file, 'w')
for contigIdx in range(len(contig_labels_dict)):
    f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
f.close()


res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                                weights=weis[0], initial_membership = initial_list,
                                     node_sizes=length_weight)


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res1 , is_membership_fixed=is_membership_fixed,
                                                  n_iterations=-1)


part = list(res1)


res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                                weights=weis[1], initial_membership = initial_list,
                                     node_sizes=length_weight)
optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res2 , is_membership_fixed=is_membership_fixed,
                                        n_iterations=-1)


part = list(res2)

###chuanlian


res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                     weights=weis[1], initial_membership = initial_list,
                                     node_sizes=length_weight)
optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res2, is_membership_fixed=is_membership_fixed,
                                        n_iterations=-1)


res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                     initial_membership =res2.membership,weights=weis[0],
                                     node_sizes=length_weight)


optimiser = leidenalg.Optimiser()

diff = 1
while diff > 0:
    diff = optimiser.optimise_partition(res1 , is_membership_fixed=is_membership_fixed,
                                        n_iterations=-1)


part = list(res1)

max_edges = 50
p1 = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10)
p2 = fit_hnsw_index(logger, norm_embeddings2, ef=max_edges * 10)



def cluster_hnswlib_leiden_markerseed_RBERVertexPartition_nodesize_partgraph_chuanlian(logger, length_weight, seed_file, max_edges, norm_embeddings, p, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None, resolution_parameter=1.0,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    graphs = []
    weis = []
    time_start = time.time()
    ann_neighbor_indices, ann_distances = p1.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()


    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff
    # each_node_index = np.arange(len(wei))
    #
    # #save at least one edge for each node
    # each_node_index = [i % max_edges == 0 for i in each_node_index]
    #
    # save_index_all = [each_node_index[i] or save_index[i] for i in range(len(save_index))]
    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    logger.info("bandwidth:\t" +str(bandwidth))

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei)


    ### coverage graph
    time_start = time.time()
    ann_neighbor_indices, ann_distances = p2.knn_query(norm_embeddings2, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources = np.repeat(np.arange(len(norm_embeddings2)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()


    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    logger.info("bandwidth:\t" +str(bandwidth))

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings2)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei)

    res2 = leidenalg.RBERVertexPartition(graphs[1],
                                             weights=weis[1], initial_membership = initial_list,
                                             node_sizes=length_weight,resolution_parameter=resolution_parameter)
    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition(res2, is_membership_fixed=is_membership_fixed,
                                            n_iterations=-1)


    res1 = leidenalg.RBERVertexPartition(graphs[0],
                                             initial_membership =res2.membership, weights=weis[0],
                                             node_sizes=length_weight,resolution_parameter=resolution_parameter)


    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition(res1, is_membership_fixed=is_membership_fixed,
                                            n_iterations=-1)


    part = list(res1)

    return part, bandwidth


#
####cami_gt
output_path='/home/wzy/data/cami_gt/output/COMEBin_whole_addcovloss_batchsize_512_test_naug_mean_embeddings/'
seed_file='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
emb_file=output_path+'embeddings.tsv'
contig_file='/home/wzy/data/cami_gt/input/gsa_f1k_atgc.fasta'


##cami_mousegut
output_path='/home/wzy/data/cami_mousegut/output/COMEBin_whole_addcovloss_batchsize_512/'
seed_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr19_09-30-42_ZzStudio-7048-4x1080/input/anonymous_gsa_pooled_f1k.fasta.bacar_marker.2quarter_lencutoff_1001.seed'
emb_file=output_path+'embeddings.tsv'
contig_file='/home/wzy/tools/mycode_gitlab_comebinner/COMEBin/CLmodel/runs/runs_trainingdatasets/Apr19_09-30-42_ZzStudio-7048-4x1080/input/anonymous_gsa_pooled_f1k.fasta'


emb_file2=output_path+'naug_mean_embeddings.tsv'


output_path = output_path + '/cluster_res/'

embHeader = pd.read_csv(emb_file, sep='\t', nrows=1)
embMat = pd.read_csv(emb_file, sep='\t', usecols=range(1, embHeader.shape[1])).values
namelist = pd.read_csv(emb_file, sep='\t', usecols=range(1)).values[:, 0]

lengths = get_length(contig_file)
length_weight = []
for seq_id in namelist:
    length_weight.append(lengths[seq_id])

norm_embeddings = normalize(embMat)

embHeader2 = pd.read_csv(emb_file2, sep='\t', nrows=1)
embMat2 = pd.read_csv(emb_file2, sep='\t', usecols=range(1, embHeader2.shape[1])).values
norm_embeddings2 = normalize(embMat2)

norm_embeddings = np.vstack((norm_embeddings, norm_embeddings2))

length_weight = length_weight * 2


max_edges = 50
p = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10)

def cluster_hnswlib_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_addmean(logger, length_weight, seed_file, max_edges, norm_embeddings, p, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    print(seed_bacar_marker_idx)

    initial_list = list(np.arange(int(len(namelist))))
    fixed_list = [i in seed_bacar_marker_idx for i in initial_list]

    initial_list = initial_list * 2
    is_membership_fixed = fixed_list*2

    time_start = time.time()
    ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()

    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    logger.info("bandwidth:\t" +str(bandwidth))

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    res = leidenalg.SurpriseVertexPartition(g,
                                            weights=wei, initial_membership = initial_list,
                                            node_sizes=length_weight)

    optimiser = leidenalg.Optimiser()
    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed,n_iterations=-1)


    part = list(res)

    return part, bandwidth


namelist = list(namelist)
namelist = namelist + [key+'_mean' for key in namelist]

# output_file = output_path+'leiden_markerseed_SurpriseVertexPartition_nodesize_nokeep1_hnsw_l2mode_bandwidth_0.2_res_maxedges50_partgraph_50.tsv.addmean.multiplex.tsv'
output_file = output_path+'leiden_markerseed_SurpriseVertexPartition_nodesize_nokeep1_hnsw_l2mode_bandwidth_0.2_res_maxedges50_partgraph_50.tsv.onlyaugmean.tsv'


contig_labels_dict ={}
# dict of communities
numnode = 0
rang = []
for ci in range(len(part)):
    rang.append(ci)
    numnode = numnode+len(part[ci])
    for id in part[ci]:
        contig_labels_dict[namelist[id]] = 'group'+str(ci)

print(output_file)
f = open(output_file, 'w')
for contigIdx in range(len(contig_labels_dict)):
    f.write(namelist[contigIdx] + "\t" + str(contig_labels_dict[namelist[contigIdx]]) + "\n")
f.close()



def cluster_hnswlib_leiden_markerseed_SurpriseVertexPartition_nodesize_partgraph_multiplex(logger, length_weight, seed_file, max_edges, norm_embeddings, p, namelist, num_process=40, percentile=5, bandwidth=None, lmode='l1', initial_membership=None, resolution_parameter=1.0,add_lengthweight=False, resultfile=None,partgraph_ratio=50 ):
    seed_bacar_marker_idx = gen_seed_idx(seed_file, contig_id_list=namelist)
    initial_list = list(np.arange(len(namelist)))
    is_membership_fixed = [i in seed_bacar_marker_idx for i in initial_list]

    graphs = []
    weis = []
    time_start = time.time()
    ann_neighbor_indices, ann_distances = p1.knn_query(norm_embeddings, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()


    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    logger.info("bandwidth:\t" +str(bandwidth))

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei)


    ### coverage graph
    time_start = time.time()
    ann_neighbor_indices, ann_distances = p2.knn_query(norm_embeddings2, max_edges+1)
    #ann_distances是l2 distance的平方
    time_end = time.time()
    logger.info('knn query time cost:\t' +str(time_end - time_start) + "s")
    sources = np.repeat(np.arange(len(norm_embeddings2)), max_edges)
    targets_indices = ann_neighbor_indices[:,1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:,1:]
    wei = wei.flatten()


    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    if lmode == 'l1':
        wei = np.sqrt(wei)
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    if lmode == 'l2':
        if not bandwidth:
            bandwidth = np.percentile(wei, percentile)
        wei = np.exp(-wei / bandwidth)

    logger.info("bandwidth:\t" +str(bandwidth))

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings2)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    graphs.append(g)
    weis.append(wei)

    res2 = leidenalg.SurpriseVertexPartition(graphs[1],
                                         weights=weis[1], initial_membership = initial_list,
                                         node_sizes=length_weight)

    res1 = leidenalg.SurpriseVertexPartition(graphs[0],
                                             initial_membership = initial_list, weights=weis[0],
                                         node_sizes=length_weight)


    optimiser = leidenalg.Optimiser()

    diff = 1
    while diff > 0:
        diff = optimiser.optimise_partition_multiplex([res1,res2], is_membership_fixed=is_membership_fixed,
                                            n_iterations=-1)


    part = list(res1)

    return part, bandwidth