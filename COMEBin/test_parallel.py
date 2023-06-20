def run_leiden(output_file, namelist,
               ann_neighbor_indices, ann_distances,length_weight, max_edges, norm_embeddings, percentile=5, bandwidth=None, lmode='l2', initial_list=None,is_membership_fixed=None, resolution_parameter=1.0,partgraph_ratio=50 ):
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

    index = sources > targets
    sources = sources[index]
    targets = targets[index]
    wei = wei[index]
    vcount = len(norm_embeddings)
    edgelist = list(zip(sources, targets))
    g = Graph(vcount, edgelist)

    res = leidenalg.RBERVertexPartition(g,
                                        weights=wei, initial_membership = initial_list,
                                        resolution_parameter = resolution_parameter,node_sizes=length_weight)

    optimiser = leidenalg.Optimiser()
    optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed,n_iterations=-1)

    # diff = 1
    # while diff > 0:
    #     diff = optimiser.optimise_partition(res, is_membership_fixed=is_membership_fixed,n_iterations=-1)

    part = list(res)


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





