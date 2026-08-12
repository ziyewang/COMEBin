# modified from https://github.com/sthalles/SimCLR
from contextlib import nullcontext
import logging
import os
from pathlib import Path
import time

import torch
import torch.nn.functional as F
from torch.amp import GradScaler, autocast
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm
from embedding_bundle import invalidate_embedding_manifest, publish_embedding_bundle
from utils import _process_memory_status, save_config_file, accuracy, save_checkpoint
import numpy as np


class SimCLR(object):

    def __init__(self, *args, **kwargs):
        """
        Initialize the SimCLR model and related components.

        :param kwargs: Keyword arguments including 'args', 'model', 'optimizer', 'scheduler'.
        """
        self.args = kwargs['args']
        self.logger = kwargs.get('logger', logging.getLogger(__name__))
        self.use_amp = self.args.fp16_precision
        parameter_count = sum(parameter.numel() for parameter in kwargs['model'].parameters())
        self.logger.info("NN setup: move model to %s; parameters=%d.", self.args.device, parameter_count)
        self.model = kwargs['model'].to(self.args.device)
        if self.args.device.type == 'cuda':
            self.logger.info("NN setup: model is on GPU; allocated=%d MiB, reserved=%d MiB.",
                             torch.cuda.memory_allocated() // 1024**2,
                             torch.cuda.memory_reserved() // 1024**2)
        self.optimizer = kwargs['optimizer']
        self.scheduler = kwargs['scheduler']
        self.logger.info("NN setup: create TensorBoard writer in %s.", self.args.output_path)
        self.writer = SummaryWriter(log_dir=self.args.output_path)
        logging.basicConfig(filename=os.path.join(self.args.output_path, 'training.log'), level=logging.DEBUG)
        self.criterion = torch.nn.CrossEntropyLoss().to(self.args.device)
        self.logger.info("NN setup finished; loss function and training state are ready.")

    def _autocast_context(self):
        if self.use_amp:
            return autocast(device_type='cuda', dtype=torch.float16)
        return nullcontext()

    def _gpu_memory_status(self):
        if self.args.device.type != 'cuda':
            return None
        return {
            'allocated_mib': torch.cuda.memory_allocated(self.args.device) // 1024**2,
            'reserved_mib': torch.cuda.memory_reserved(self.args.device) // 1024**2,
            'peak_allocated_mib': torch.cuda.max_memory_allocated(self.args.device) // 1024**2,
        }

    def _log_training_start(self, mode, train_loader, epochs):
        self.logger.info(
            "NN training preflight: mode=%s, epochs=%d, samples=%d, batches_per_epoch=%d, "
            "batch_size=%d, device=%s, fp16=%s, process_memory=%s, gpu_memory=%s.",
            mode, epochs, len(train_loader.dataset), len(train_loader), train_loader.batch_size,
            self.args.device, self.use_amp, _process_memory_status(), self._gpu_memory_status())

    def _log_epoch_progress(self, epoch_counter, total_epochs, batch_count, loss, top1, epoch_start):
        self.logger.info(
            "NN training: epoch=%d/%d, batches=%d, loss=%.8g, top1=%.4f, elapsed=%.3fs, "
            "process_memory=%s, gpu_memory=%s.",
            epoch_counter + 1, total_epochs, batch_count, loss.detach().item(), top1[0].item(),
            time.monotonic() - epoch_start, _process_memory_status(), self._gpu_memory_status())

    def info_nce_loss(self, features):
        """
        Calculate the InfoNCE loss for SimCLR.

        :param features: Input features.
        :return: Logits and labels for the loss.
        """
        sample_ids = torch.arange(self.args.batch_size, device=features.device).repeat(self.args.n_views)
        positive_mask = sample_ids.unsqueeze(0).eq(sample_ids.unsqueeze(1))

        features = F.normalize(features, dim=1)

        similarity_matrix = torch.matmul(features, features.T)

        # discard the main diagonal from both: labels and similarities matrix
        diagonal_mask = torch.eye(positive_mask.shape[0], dtype=torch.bool, device=features.device)
        positive_mask = positive_mask[~diagonal_mask].view(positive_mask.shape[0], -1)
        similarity_matrix = similarity_matrix[~diagonal_mask].view(similarity_matrix.shape[0], -1)

        # select and combine multiple positives
        positives = similarity_matrix[positive_mask].view(-1, 1)
        # select only the negatives the negatives
        negatives = similarity_matrix[~positive_mask].view(similarity_matrix.shape[0], -1)
        negatives = negatives[:, None].expand(-1, self.args.n_views - 1, -1).flatten(0, 1)

        logits = torch.cat([positives, negatives], dim=1)
        labels = torch.zeros(logits.shape[0], dtype=torch.long, device=features.device)

        logits = logits / self.args.temperature
        return logits, labels

    def _publish_embeddings(self, data, namelist, contig_lengths, kmer_len):
        output_path = Path(self.args.output_path)
        n_contigs = len(data[0])
        if len(namelist) != n_contigs or len(contig_lengths) != n_contigs:
            raise ValueError(
                f'Embedding publication axes are not aligned: features={n_contigs}, '
                f'identifiers={len(namelist)}, lengths={len(contig_lengths)}.')
        invalidate_embedding_manifest(output_path)
        embedding_temp = output_path / f'.embeddings.{os.getpid()}.tmp.npy'
        coverage_temp = output_path / f'.covembeddings.{os.getpid()}.tmp.npy'
        embedding_temp.unlink(missing_ok=True)
        coverage_temp.unlink(missing_ok=True)
        embedding_target = None
        coverage_target = None
        self.logger.info(
            'Final inference: start one traversal; contigs=%d, batch_size=1024, '
            'embedding_dimensions=%d, coverage_dimensions=%d.',
            n_contigs, self.args.out_dim, self.args.out_dim_forcov)
        try:
            embedding_target = np.lib.format.open_memmap(
                embedding_temp, mode='w+', dtype='<f4',
                shape=(n_contigs, self.args.out_dim))
            coverage_target = np.lib.format.open_memmap(
                coverage_temp, mode='w+', dtype='<f4',
                shape=(n_contigs, self.args.out_dim_forcov))
            with torch.no_grad():
                self.model.eval()
                for start in range(0, n_contigs, 1024):
                    end = min(start + 1024, n_contigs)
                    kmer_batch = data[0][start:end, -kmer_len:].to(self.args.device)
                    coverage_batch = data[0][start:end, :-kmer_len].to(self.args.device)
                    model_outputs = self.model(kmer_batch, coverage_batch)
                    if not isinstance(model_outputs, tuple) or len(model_outputs) < 2:
                        raise ValueError(
                            'The combined model did not return combined and coverage embeddings.')
                    embeddings = model_outputs[0].detach().to('cpu').numpy()
                    coverage_embeddings = model_outputs[1].detach().to('cpu').numpy()
                    expected_rows = end - start
                    if (embeddings.dtype != np.float32 or
                            embeddings.shape != (expected_rows, self.args.out_dim) or
                            not np.isfinite(embeddings).all()):
                        raise ValueError(
                            f'Invalid combined embeddings for rows {start}:{end}: '
                            f'dtype={embeddings.dtype}, shape={embeddings.shape}.')
                    if (coverage_embeddings.dtype != np.float32 or
                            coverage_embeddings.shape != (expected_rows, self.args.out_dim_forcov) or
                            not np.isfinite(coverage_embeddings).all()):
                        raise ValueError(
                            f'Invalid coverage embeddings for rows {start}:{end}: '
                            f'dtype={coverage_embeddings.dtype}, '
                            f'shape={coverage_embeddings.shape}.')
                    embedding_target[start:end] = embeddings
                    coverage_target[start:end] = coverage_embeddings
                    if start == 0:
                        self.logger.info(
                            'Final inference: first batch published to temporary targets; '
                            'rows=%d, process_memory=%s, gpu_memory=%s.',
                            expected_rows, _process_memory_status(), self._gpu_memory_status())
            embedding_target.flush()
            coverage_target.flush()
        except BaseException:
            if embedding_target is not None:
                del embedding_target
            if coverage_target is not None:
                del coverage_target
            embedding_temp.unlink(missing_ok=True)
            coverage_temp.unlink(missing_ok=True)
            raise
        else:
            del embedding_target, coverage_target

        publish_embedding_bundle(
            output_path, embedding_temp, coverage_temp, namelist,
            contig_lengths, self.args.out_dim, self.args.out_dim_forcov)
        self.logger.info(
            'Final inference: published complete binary embedding bundle; '
            'contigs=%d.', n_contigs)

    def train(self, train_loader, data, namelist, contig_lengths):
        """
        Train the combined k-mer and coverage model.

        :param train_loader: Data loader for training.
        :param data: Input data.
        :param namelist: List of sequence names.
        """
        scaler = GradScaler(self.args.device.type, enabled=self.use_amp)

        # save config file
        save_config_file(self.args.output_path, self.args)

        earlystop_epoch=0
        logging.info(f"Start SimCLR training for {self.args.epochs} epochs.")
        self._log_training_start('combined-model', train_loader, self.args.epochs)
        self.logger.info("NN training: enter combined-model loop and request the first DataLoader batch.")
        kmer_len = 136
        logging.info('kmer_len:\t' + str(kmer_len) + '\n')

        for epoch_counter in range(self.args.epochs):
            epoch_start = time.monotonic()
            batch_count = 0

            for batch_index, contig_features in enumerate(tqdm(train_loader)):
                batch_count += 1
                if epoch_counter == 0 and batch_index == 0:
                    self.logger.info("NN training: first combined-model batch received; view_shapes=%s.",
                                     [tuple(view.shape) for view in contig_features])
                contig_features = torch.cat(contig_features, dim=0)

                contig_features = contig_features.to(self.args.device)
                if epoch_counter == 0 and batch_index == 0:
                    self.logger.info("NN training: first combined-model batch moved to %s; shape=%s.",
                                     self.args.device, tuple(contig_features.shape))
                with self._autocast_context():
                    features, _ = self.model(contig_features[:, -kmer_len:], contig_features[:, :-kmer_len])
                    logits, labels = self.info_nce_loss(features)
                    loss = self.criterion(logits, labels)
                if epoch_counter == 0 and batch_index == 0:
                    self.logger.info("NN training: first combined-model forward pass finished; loss=%s.",
                                     loss.detach().item())

                self.optimizer.zero_grad()

                scaler.scale(loss).backward()

                scaler.step(self.optimizer)
                scaler.update()
                if epoch_counter == 0 and batch_index == 0:
                    self.logger.info("NN training: first combined-model backward pass and optimizer step finished.")

            if batch_count == 0:
                raise RuntimeError(
                    f"Epoch {epoch_counter + 1} received no DataLoader batches; "
                    f"samples={len(train_loader.dataset)}, batch_size={train_loader.batch_size}, "
                    f"drop_last={train_loader.drop_last}.")
            top1, top5 = accuracy(logits, labels, topk=(1, min(5, logits.shape[1])))
            self._log_epoch_progress(
                epoch_counter, self.args.epochs, batch_count, loss, top1, epoch_start)

            if not self.args.notuse_scheduler:
                # warmup for the first 10 epochs
                if epoch_counter >= 10:
                    self.scheduler.step()
            logging.debug(f"Epoch: {epoch_counter}\tLoss: {loss}\tTop1 accuracy: {top1[0]}")

            if self.args.earlystop:
                if epoch_counter >= 10 and top1[0] > 99.0:
                    earlystop_epoch +=1
                else:
                    earlystop_epoch = 0

                if earlystop_epoch >=3:
                    break

        logging.info("Training has finished.")
        # save model checkpoints
        checkpoint_name = 'checkpoint_{:04d}.pth.tar'.format(self.args.epochs)
        save_checkpoint({
            'epoch': self.args.epochs,
            # 'arch': self.args.arch,
            'state_dict': self.model.state_dict(),
            'optimizer': self.optimizer.state_dict(),
        }, is_best=False, filename=os.path.join(self.args.output_path, checkpoint_name))
        logging.info(f"Model checkpoint and metadata has been saved at {self.args.output_path}.")

        # ckpt = torch.load('/checkpoint_0200.pth.tar')

        self._publish_embeddings(data, namelist, contig_lengths, kmer_len)
