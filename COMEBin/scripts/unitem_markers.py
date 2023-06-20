# from https://github.com/dparks1134/UniteM/blob/master/unitem/markers.py
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.markerSets import MarkerSetParser, BinMarkerSets

from scripts.unitem_defaults import *

"""
CHECKM_BAC_MS = 'bacteria.ms'
CHECKM_AR_MS = 'archaea.ms'
"""


class Markers():
    """Read marker information."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        parser = MarkerSetParser()
        bin_marker_sets = parser.parseTaxonomicMarkerSetFile(CHECKM_BAC_MS)
        self.bac_ms = bin_marker_sets.mostSpecificMarkerSet()

        bin_marker_sets = parser.parseTaxonomicMarkerSetFile(CHECKM_AR_MS)
        self.ar_ms = bin_marker_sets.mostSpecificMarkerSet()

        self.bac_markers_on_contigs = None
        self.ar_markers_on_contigs = None

    # update for metabinner
    def read_table(self, marker_gene_table):
        """Read table indicating marker genes for each bin."""

        markers = defaultdict(lambda: defaultdict(list))
        if not os.path.exists(marker_gene_table):
            # did not identify any markers
            return markers

        with open(marker_gene_table) as f:
            f.readline()

            for line in f:
                line_split = line.strip().split('\t')

                # bin_id = line_split[0]
                marker_id = line_split[1]
                gene_id = line_split[2]
                if '&&' in gene_id:
                    # indicates a marker gene is identified
                    # in adjacent genes and is likely an
                    # assembly or gene calling error
                    gene_ids = gene_id.split('&&')
                    gene_id = gene_ids[0]
                scaffold_id = gene_id[0:gene_id.rfind('_')]

                if marker_id not in markers:
                    markers[marker_id] = [scaffold_id]
                else:
                    markers[marker_id].append(scaffold_id)

        return markers

    def evaluate_bac(self, gene_table, individual_markers=False):
        """Evaluate completeness and contamination of genome using bacterial marker sets."""

        comp, cont = self.bac_ms.genomeCheck(gene_table, individual_markers)
        return comp, cont

    def evaluate_ar(self, gene_table, individual_markers=False):
        """Evaluate completeness and contamination of genome using archaeal marker sets."""

        comp, cont = self.ar_ms.genomeCheck(gene_table, individual_markers)
        return comp, cont

    def evaluate(self, bac_gene_table, ar_gene_table):
        """Evaluate completeness and contamination of genome using best domain-level marker sets."""

        bac_comp, bac_cont = self.evaluate_bac(bac_gene_table, True)
        ar_comp, ar_cont = self.evaluate_ar(ar_gene_table, True)

        if bac_comp + bac_cont > ar_comp + ar_cont:
            # select domain set with the larget number of identified markers
            # including those present multiple times
            bac_comp, bac_cont = self.evaluate_bac(bac_gene_table, False)
            return 'Bacteria', bac_comp, bac_cont

        ar_comp, ar_cont = self.evaluate_ar(ar_gene_table, False)
        return 'Archaea', ar_comp, ar_cont

    def bin_quality(self, bin):
        """Estimate quality of bin."""

        # create gene tables for bin
        bac_gene_table = defaultdict(list)
        ar_gene_table = defaultdict(list)
        for cid in bin:
            for marker_id in self.bac_markers_on_contigs[cid]:
                bac_gene_table[marker_id].append(cid)

            for marker_id in self.ar_markers_on_contigs[cid]:
                ar_gene_table[marker_id].append(cid)

        domain, comp, cont = self.evaluate(bac_gene_table,
                                           ar_gene_table)

        return domain, comp, cont

    # changed for metabinner
    def marker_gene_tables(self, bac_mg_table, ar_mg_table):
        """Get marker genes for bins across all binning methods."""

        markers = Markers()
        # binning_methods_dir = os.path.join(profile_dir, BINNING_METHOD_DIR)

        # bac_mg_table = os.path.join(binning_methods_dir, bm, CHECKM_BAC_DIR, MARKER_GENE_TABLE)
        # bac_mg_table ='/home1/wangzy/binning_data/sim40_20/output/metabinner_res/unitem/my_profile/binning_methods/X_cov_notrans1/checkm_bac/marker_gene_table.tsv'
        bac_gene_tables = markers.read_table(bac_mg_table)

        # ar_mg_table = os.path.join(binning_methods_dir, bm, CHECKM_AR_DIR, MARKER_GENE_TABLE)
        # ar_mg_table ='/home1/wangzy/binning_data/sim40_20/output/metabinner_res/unitem/my_profile/binning_methods/X_cov_notrans1/checkm_ar/marker_gene_table.tsv'
        ar_gene_tables = markers.read_table(ar_mg_table)

        gene_tables = (bac_gene_tables, ar_gene_tables)

        self._markers_on_contigs(gene_tables)

        return gene_tables

    # update for metabinner
    def _markers_on_contigs(self, gene_tables):
        """Get markers on each contig."""

        self.bac_markers_on_contigs = defaultdict(list)
        self.ar_markers_on_contigs = defaultdict(list)
        processed_contigs = set()
        (bac_gene_tables, ar_gene_tables) = gene_tables
        scaffolds_in_binning_method = set()
        for marker_id, scaffold_ids in bac_gene_tables.items():
            for scaffold_id in scaffold_ids:
                if scaffold_id in processed_contigs:
                    continue
                self.bac_markers_on_contigs[scaffold_id].append(marker_id)
            scaffolds_in_binning_method.update(scaffold_ids)

        for marker_id, scaffold_ids in ar_gene_tables.items():
            for scaffold_id in scaffold_ids:
                if scaffold_id in processed_contigs:
                    continue
                self.ar_markers_on_contigs[scaffold_id].append(marker_id)
            scaffolds_in_binning_method.update(scaffold_ids)

        processed_contigs.update(scaffolds_in_binning_method)

    def create_gene_table(self, contig_ids):
        """Create gene tables for a set of contigs."""

        bac_table = defaultdict(list)
        ar_table = defaultdict(list)
        for cid in contig_ids:
            for mid in self.bac_markers_on_contigs[cid]:
                bac_table[mid].append(cid)

            for mid in self.ar_markers_on_contigs[cid]:
                ar_table[mid].append(cid)

        return bac_table, ar_table
