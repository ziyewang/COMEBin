import torch
import torch.backends.cudnn as cudnn
import torch.nn as nn
import numpy as np
import os
import sys

from get_augfeature import get_ContrastiveLearningDataset
from utils import get_length

from models.mlp import EmbeddingNet
from simclr import SimCLR

torch.manual_seed(1)
torch.cuda.manual_seed_all(1)
# torch.cuda.manual_seed(0)
torch.backends.cudnn.deterministic = True


def train_CLmodel(logger, args):
    if torch.cuda.is_available():
        args.device = torch.device('cuda')
        cudnn.deterministic = True
        cudnn.benchmark = True
    else:
        args.device = torch.device('cpu')

    logger.info("Generate features for the data.")

    if os.path.exists(args.output_path+'/embeddings.tsv'):
        logger.info("The embeddings file has been generated before, please check the output commands.")
        sys.exit()

    dataset, namelist = get_ContrastiveLearningDataset(args.data, args.n_views,
                                                       args.kmer_model_path, args.device, args.nokmer, args.cov_meannormalize,
                                                       args.cov_minmaxnormalize, args.cov_standardization,args.addvars,args.vars_sqrt, args.kmer_l2_normalize, args.kmerMetric_notl2normalize)

    contig_file = args.data + '/aug0/sequences_aug0.fasta'
    lengths = get_length(contig_file)
    length_weight = []
    for seq_id in namelist:
        length_weight.append(lengths[seq_id])


    if args.notuseaug0:
        args.n_views = args.n_views - 1
        train_dataset = torch.utils.data.TensorDataset(*[dataset[i+1][np.array(length_weight) >= args.contig_len]
                                                         for i in range(args.n_views)])
    else:
        # train_dataset = torch.utils.data.TensorDataset(*dataset)
        train_dataset = torch.utils.data.TensorDataset(*[dataset[i][np.array(length_weight) >= args.contig_len]
                                                         for i in range(args.n_views)])

    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True,
        num_workers=args.workers, pin_memory=True, drop_last=True)

    # Set embedder model.
    if not args.add_model_for_coverage:
        ps = [args.dropout_value] * (args.n_layer - 1)
        actn = nn.LeakyReLU()

        emb_szs_list = [args.emb_szs] * args.n_layer

        model = EmbeddingNet(
            in_sz=len(dataset[0][0]),
            out_sz=args.out_dim,
            emb_szs=emb_szs_list,
            ps=ps,
            use_bn=True,
            actn=actn,
        )

        optimizer = torch.optim.AdamW(model.parameters(), args.lr, weight_decay=args.weight_decay)

        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=0,
                                                               last_epoch=-1)

        simclr = SimCLR(model=model, optimizer=optimizer, scheduler=scheduler, args=args)
        simclr.train(train_loader, dataset, namelist)

    else:
        if args.pretrain_kmer_model_path !='no':
            config_file= os.path.dirname(args.pretrain_kmer_model_path)+'/kmerMetric_config.yaml'

            from ruamel.yaml import YAML
            from pathlib import Path

            yaml = YAML(typ='safe')

            cnf = yaml.load(Path(config_file))

            ps = [cnf['dropout_value']]*(len(cnf['emb_szs'])-1)
            actn= nn.LeakyReLU()

            insize = 136

            kmerMetric_model = EmbeddingNet(
                in_sz=insize,
                out_sz=cnf['embedding_dim'],
                emb_szs=cnf['emb_szs'],
                ps=ps,
                use_bn=True,
                actn=actn,
            )

            kmerMetric_model = kmerMetric_model.to(args.device)

            if not args.not_load_kmermetric_state:
                kmerMetric_model.load_state_dict(torch.load(args.pretrain_kmer_model_path, map_location=args.device))

            cov_dim = len(dataset[0][0]) - insize
            input_size = args.out_dim_forcov + 128
            print('cov_dim:\t' + str(cov_dim) + '\n')

            emb_szs_list = [args.emb_szs_forcov] * args.n_layer_forcov

            cov_model = EmbeddingNet(
                in_sz=cov_dim,
                out_sz=args.out_dim_forcov,
                emb_szs=emb_szs_list,
                ps=[args.dropout_value] * (args.n_layer_forcov - 1),
                use_bn=True,
                actn=nn.LeakyReLU(),
            )


            from models.mlp2 import EmbeddingNet as EmbeddingNet2

            emb_szs_list = [args.emb_szs] * args.n_layer

            model = EmbeddingNet2(
                in_sz=input_size,
                out_sz=args.out_dim,
                emb_szs=emb_szs_list,
                ps=[args.dropout_value] * (args.n_layer - 1),
                use_bn=True,
                actn=nn.LeakyReLU(),
                cov_model=cov_model,
                pretrained_model=kmerMetric_model,
                covmodel_notl2normalize=args.covmodel_notl2normalize,
            )

            if args.finetunepretrainmodel:
                pretrained_model_paraname = []
                for name, param in model.pretrained_model.named_parameters():
                    # print(name)
                    pretrained_model_paraname.append('pretrained_model.'+name)
                print(pretrained_model_paraname)

                params_lx = [param for name, param in model.named_parameters()
                             if name not in pretrained_model_paraname]
                # print(params_lx)
                # for name, param in model.named_parameters():
                #     print(name)
                # # print(model.named_parameters())
                # print(model.pretrained_model)
                # os._exit()
                optimizer = torch.optim.AdamW([{'params': params_lx},
                                               {'params': model.pretrained_model.parameters(),
                                                'lr': args.lr * args.finetunelr_ratio}],
                                              lr=args.lr, weight_decay=args.weight_decay)


            else:
                optimizer = torch.optim.AdamW(model.parameters(), args.lr, weight_decay=args.weight_decay)

            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=0,
                                                                   last_epoch=-1)

            simclr = SimCLR(model=model, optimizer=optimizer, scheduler=scheduler, args=args)
            simclr.train_addpretrain(train_loader, dataset, namelist)



        else:
            if args.kmer_model_path == 'empty':
                cov_dim = len(dataset[0][0]) - 136
                input_size = args.out_dim_forcov + 136
            else:
                cov_dim = len(dataset[0][0]) - 128
                input_size = args.out_dim_forcov + 128
                print('cov_dim:\t' + str(cov_dim) + '\n')


            emb_szs_list = [args.emb_szs_forcov] * args.n_layer_forcov

            cov_model = EmbeddingNet(
                in_sz=cov_dim,
                out_sz=args.out_dim_forcov,
                emb_szs=emb_szs_list,
                ps=[args.dropout_value] * (args.n_layer_forcov - 1),
                use_bn=True,
                actn=nn.LeakyReLU(),
            )

            if args.pretrain_coveragemodel:
                print("pretrain_coveragemodel!")
                print(cov_model.state_dict())
                print(len(dataset[0][0][128:]))
                optimizer = torch.optim.AdamW(cov_model.parameters(), args.lr, weight_decay=args.weight_decay)

                scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.covmodelepochs, eta_min=0,
                                                                       last_epoch=-1)

                simclr = SimCLR(model=cov_model, optimizer=optimizer, scheduler=scheduler, args=args)
                simclr.covmodeltrain(train_loader)
                print(cov_model.state_dict())
                print(len(dataset[0][0][128:]))

            from models.mlp2 import EmbeddingNet as EmbeddingNet2


            emb_szs_list = [args.emb_szs] * args.n_layer

            model = EmbeddingNet2(
                in_sz=input_size,
                out_sz=args.out_dim,
                emb_szs=emb_szs_list,
                ps=[args.dropout_value] * (args.n_layer - 1),
                use_bn=True,
                actn=nn.LeakyReLU(),
                cov_model=cov_model,
                covmodel_notl2normalize=args.covmodel_notl2normalize,
            )
            optimizer = torch.optim.AdamW(model.parameters(), args.lr, weight_decay=args.weight_decay)

            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs, eta_min=0,
                                                                   last_epoch=-1)

            simclr = SimCLR(model=model, optimizer=optimizer, scheduler=scheduler, args=args)
            simclr.train_addpretrain(train_loader, dataset, namelist)

    logger.info("Finish training.")
