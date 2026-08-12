import os
import re
import sys
from pathlib import Path

import torch
import torch.nn as nn
import numpy as np

from get_augfeature import get_ContrastiveLearningDataset
from utils import _process_memory_status, configure_reproducibility, seed_data_loader_worker

from models.mlp import EmbeddingNet
from models.mlp2 import CombinedEmbeddingNet
from simclr import SimCLR


DATALOADER_WORKERS = 4


def _resolve_training_device(requested):
    requested = requested.strip().lower()
    if requested == 'cpu':
        return torch.device('cpu')
    if requested != 'cuda' and re.fullmatch(r'cuda:\d+', requested) is None:
        raise ValueError(f"Unsupported device {requested!r}; expected cuda, cuda:N, or cpu.")
    cuda_environment = (
        f"PyTorch CUDA={torch.version.cuda}, "
        f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
    if torch.version.cuda is None:
        raise RuntimeError(
            "CUDA training was requested, but the installed PyTorch is CPU-only; "
            f"{cuda_environment}. Use --device cpu only if CPU training is intended.")
    try:
        cuda_available = torch.cuda.is_available()
        visible_devices = torch.cuda.device_count()
    except Exception as error:
        raise RuntimeError(
            f"CUDA discovery failed for requested device {requested!r}: {error}; "
            f"{cuda_environment}.") from error
    if not cuda_available:
        raise RuntimeError(
            "CUDA training was requested, but CUDA is unavailable; "
            f"{cuda_environment}, visible_devices={visible_devices}.")

    requested_device = torch.device(requested)
    try:
        device_index = torch.cuda.current_device() if requested_device.index is None else requested_device.index
    except Exception as error:
        raise RuntimeError(
            f"CUDA could not select the current device for {requested!r}: {error}; "
            f"{cuda_environment}, visible_devices={visible_devices}.") from error
    if device_index >= visible_devices:
        raise RuntimeError(
            f"Requested cuda:{device_index}, but only {visible_devices} CUDA devices are visible; "
            f"{cuda_environment}.")

    device = torch.device('cuda', device_index)
    try:
        torch.cuda.set_device(device)
        probe = torch.ones(1, device=device)
        probe.add_(1).sum().item()
        torch.cuda.synchronize(device)
    except Exception as error:
        raise RuntimeError(
            f"CUDA device {device} is visible but failed its initialization test: {error}; "
            f"{cuda_environment}, visible_devices={visible_devices}.") from error
    return device


def train_CLmodel(logger, args):
    """
    Train the Contrastive Learning model.
    """
    if args.seed is not None:
        configure_reproducibility(args.seed)
    data_path = Path(args.data).expanduser().resolve()
    output_path = Path(args.output_path).expanduser().resolve()

    if not data_path.is_dir():
        raise FileNotFoundError(f"Feature directory does not exist: {data_path}")

    if args.num_threads < 1:
        raise ValueError('num_threads must be greater than zero.')
    args.device = _resolve_training_device(args.device)
    if args.fp16_precision and args.device.type != 'cuda':
        raise ValueError('--fp16-precision requires a CUDA training device.')
    torch.set_num_threads(args.num_threads)
    logger.info(
        "Training resources: device=%s, PyTorch intra-op threads=%d, DataLoader workers=%d (fixed). KMeans, NumPy, "
        "scikit-learn, BLAS, OpenMP, and other native runtime pools are not controlled by --num_threads.",
        args.device, torch.get_num_threads(), DATALOADER_WORKERS)

    logger.info("Training diagnostics: pid=%d, Python=%s, PyTorch=%s, CUDA=%s, cuDNN=%s.",
                os.getpid(), sys.version.split()[0], torch.__version__, torch.version.cuda,
                torch.backends.cudnn.version())
    if args.seed is not None:
        logger.info("Deterministic training: seed=%d, deterministic_algorithms=%s, cuDNN_deterministic=%s, cuDNN_benchmark=%s.",
                    args.seed, torch.are_deterministic_algorithms_enabled(),
                    torch.backends.cudnn.deterministic, torch.backends.cudnn.benchmark)
    if args.device.type == 'cuda':
        device_index = args.device.index
        properties = torch.cuda.get_device_properties(device_index)
        logger.info("GPU diagnostics: logical_device=%d, name=%s, capability=%d.%d, total_memory=%d MiB, "
                    "PyTorch_CUDA=%s, CUDA_VISIBLE_DEVICES=%s.",
                    device_index, properties.name, properties.major, properties.minor,
                    properties.total_memory // 1024**2, torch.version.cuda,
                    os.environ.get('CUDA_VISIBLE_DEVICES', 'not set'))
    logger.info("Generate features: start; data=%s, views=%d, process_memory=%s.",
                data_path, args.n_views, _process_memory_status())

    if (output_path / 'embedding_manifest.json').is_file():
        from embedding_bundle import load_embedding_bundle
        try:
            load_embedding_bundle(output_path / 'embeddings.npy')
        except (OSError, ValueError) as error:
            logger.warning(
                "Existing binary embedding bundle is not reusable and will be regenerated: %s", error)
        else:
            logger.info("The binary embedding bundle has been generated before, please check the output commands.")
            sys.exit()

    logger.info("Load features from %s.", data_path)
    dataset, namelist, contig_lengths, feature_metadata = get_ContrastiveLearningDataset(
        data_path, args.n_views, args.cov_meannormalize,
        args.cov_minmaxnormalize, args.cov_standardization,
        args.addvars, args.vars_sqrt, logger=logger)
    if args.temperature is None:
        assembly_n50 = feature_metadata.get('assembly_n50')
        if isinstance(assembly_n50, bool) or not isinstance(assembly_n50, int) or assembly_n50 <= 0:
            raise ValueError(
                'Automatic temperature requires a positive integer assembly_n50 '
                'in feature_manifest.json; regenerate augmentation data or '
                'specify --temperature explicitly.')
        args.temperature = 0.07 if assembly_n50 > 10000 else 0.15
        logger.info(
            'Automatically selected temperature=%s from complete-assembly '
            'N50=%d.', args.temperature, assembly_n50)
    else:
        if not np.isfinite(args.temperature) or args.temperature <= 0:
            raise ValueError('temperature must be finite and greater than zero.')
        logger.info('Use user-specified temperature=%s.', args.temperature)
    logger.info("Generate features: finished; contigs=%d, views=%d, shapes=%s, dtypes=%s, process_memory=%s.",
                len(namelist), len(dataset), [tuple(view.shape) for view in dataset],
                [str(view.dtype) for view in dataset], _process_memory_status())

    selected_views = dataset[1:] if args.notuseaug0 else dataset
    args.n_views = len(selected_views)
    train_mask = np.asarray(contig_lengths) >= args.contig_len
    train_sample_count = int(train_mask.sum())
    if train_sample_count < 2:
        raise ValueError(
            f"NN training needs at least two contigs with length >= {args.contig_len}; "
            f"found {train_sample_count} among {len(namelist)} feature contigs.")
    if args.batch_size < 2:
        raise ValueError(f"NN training batch size must be at least two; received {args.batch_size}.")
    if args.batch_size > train_sample_count:
        logger.warning("Requested batch_size=%d exceeds the %d trainable contigs remaining after the length filter; "
                       "batch_size is reduced to %d.", args.batch_size, train_sample_count, train_sample_count)
        args.batch_size = train_sample_count

    if train_mask.all():
        training_views = selected_views
    else:
        train_mask = torch.from_numpy(train_mask)
        training_views = [view[train_mask] for view in selected_views]
    train_dataset = torch.utils.data.TensorDataset(*training_views)
    if train_sample_count % args.batch_size:
        logger.warning("DataLoader drop_last=True will discard %d of %d trainable contigs in each epoch because "
                       "batch_size=%d.", train_sample_count % args.batch_size, train_sample_count, args.batch_size)
    data_loader_generator = None
    data_loader_worker_init = None
    if args.seed is not None:
        data_loader_generator = torch.Generator()
        data_loader_generator.manual_seed(args.seed)
        data_loader_worker_init = seed_data_loader_worker
    logger.info("Build DataLoader: all_contigs=%d, trainable_contigs=%d, batch_size=%d, workers=%d, drop_last=True.",
                len(namelist), train_sample_count, args.batch_size, DATALOADER_WORKERS)
    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=args.batch_size, shuffle=True,
        num_workers=DATALOADER_WORKERS, pin_memory=args.device.type == 'cuda', drop_last=True,
        worker_init_fn=data_loader_worker_init, generator=data_loader_generator)

    logger.info("DataLoader constructed; workers start when the first batch is requested. process_memory=%s.",
                _process_memory_status())
    coverage_dimensions = dataset[0].shape[1] - 136
    if coverage_dimensions < 1:
        raise ValueError(
            f'The assembled feature tensor must contain coverage and 136 k-mer dimensions; '
            f'found {dataset[0].shape[1]} total dimensions.')
    logger.info("Initialize combined k-mer and coverage network: coverage_dimensions=%d, device=%s.",
                coverage_dimensions, args.device)
    coverage_model = EmbeddingNet(
        in_sz=coverage_dimensions, out_sz=args.out_dim_forcov,
        emb_szs=[args.emb_szs_forcov] * args.n_layer_forcov,
        ps=[args.dropout_value] * (args.n_layer_forcov - 1),
        use_bn=True, actn=nn.LeakyReLU())
    model = CombinedEmbeddingNet(
        in_sz=args.out_dim_forcov + 136, out_sz=args.out_dim,
        emb_szs=[args.emb_szs] * args.n_layer,
        ps=[args.dropout_value] * (args.n_layer - 1),
        use_bn=True, actn=nn.LeakyReLU(), cov_model=coverage_model)
    optimizer = torch.optim.AdamW(model.parameters(), args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs, eta_min=0, last_epoch=-1)
    simclr = SimCLR(model=model, optimizer=optimizer, scheduler=scheduler, args=args, logger=logger)
    simclr.train(train_loader, dataset, namelist, contig_lengths)

    logger.info("Finish training.")
