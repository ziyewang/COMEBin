"""Tests for small-data fixes in COMEBin."""
import os
import sys
import tempfile
import shutil

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))


class TestGenSeedIdx:
    """Test gen_seed_idx guard against missing/empty seed files."""

    def test_missing_file(self):
        from COMEBin.cluster import gen_seed_idx
        result = gen_seed_idx("/nonexistent/file.seed", ["contig1", "contig2"])
        assert result == [], "Should return empty list for missing file"

    def test_empty_file(self):
        from COMEBin.cluster import gen_seed_idx
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.seed')
        tmp.close()
        try:
            result = gen_seed_idx(tmp.name, ["contig1", "contig2"])
            assert result == [], "Should return empty list for empty file"
        finally:
            os.unlink(tmp.name)

    def test_no_matching_contigs(self):
        from COMEBin.cluster import gen_seed_idx
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.seed')
        tmp.write("contig3\ncontig4\n")
        tmp.close()
        try:
            result = gen_seed_idx(tmp.name, ["contig1", "contig2"])
            assert result == [], "Should return empty list when no contigs match"
        finally:
            os.unlink(tmp.name)

    def test_valid_file(self):
        from COMEBin.cluster import gen_seed_idx
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.seed')
        tmp.write("contig1\ncontig2\n")
        tmp.close()
        try:
            result = gen_seed_idx(tmp.name, ["contig1", "contig2", "contig3"])
            assert result == [0, 1], "Should return indices of matching contigs"
        finally:
            os.unlink(tmp.name)


class TestGenSeedFallback:
    """Test gen_seed graceful fallback instead of sys.exit()."""

    def test_import_sys_exit_removed(self):
        """Verify sys.exit is no longer called when tools fail."""
        import COMEBin.utils as utils
        import inspect
        source = inspect.getsource(utils.gen_seed)
        assert "sys.exit()" not in source, "gen_seed should not call sys.exit()"
        assert "candK = 0" in source, "gen_seed should fallback to candK = 0"


class TestFitHnswIndexGuard:
    """Test fit_hnsw_index guard against <2 elements."""

    def test_zero_elements(self):
        """Should return None for empty features list."""
        from COMEBin.cluster import fit_hnsw_index
        import logging
        logger = logging.getLogger('test')
        result = fit_hnsw_index(logger, [], 1)
        assert result is None, "Should return None for 0 elements"

    def test_one_element(self):
        """Should return None for single element."""
        from COMEBin.cluster import fit_hnsw_index
        import logging
        logger = logging.getLogger('test')
        result = fit_hnsw_index(logger, [[0.1]], 1)
        assert result is None, "Should return None for 1 element"


class TestClusterGuard:
    """Test cluster() guard when too few contigs after filtering."""

    def test_cluster_early_return(self):
        """Verify cluster() returns early when <2 contigs after filtering."""
        import COMEBin.cluster as cluster_mod
        import inspect
        source = inspect.getsource(cluster_mod.cluster)
        assert "if len(namelist) < 2:" in source


class TestSeedKmeansGuard:
    """Test seed_kmeans_full guard when no marker seeds found."""

    def test_empty_seeds_guard(self):
        """Verify seed_kmeans_full skips when no marker seeds."""
        import COMEBin.cluster as cluster_mod
        import inspect
        source = inspect.getsource(cluster_mod.seed_kmeans_full)
        assert "len(seed_bacar_marker_idx)" in source or "No marker seeds found" in source


class TestRunLeidenGuard:
    """Test run_leiden guard against too few edges/nodes."""

    def test_too_few_nodes(self):
        """Verify run_leiden skips when <2 nodes or 0 edges."""
        import COMEBin.cluster as cluster_mod
        import inspect
        source = inspect.getsource(cluster_mod.run_leiden)
        assert "len(sources) == 0 or vcount < 2" in source


class TestDataLoaderFix:
    """Test DataLoader drop_last fix."""

    def test_drop_last_false(self):
        """Verify drop_last is now False."""
        import COMEBin.train_CLmodel as train_mod
        import inspect
        source = inspect.getsource(train_mod.train_CLmodel)
        assert "drop_last=False" in source, "Should use drop_last=False"

    def test_empty_dataset_guard(self):
        """Verify early return when train_dataset is empty."""
        import COMEBin.train_CLmodel as train_mod
        import inspect
        source = inspect.getsource(train_mod.train_CLmodel)
        assert "len(train_dataset) == 0" in source


class TestKnnQueryCap:
    """Test kNN query k capped to available elements."""

    def test_max_k_capped(self):
        """Verify knn_query uses min(requested, available)."""
        import COMEBin.cluster as cluster_mod
        import inspect
        source = inspect.getsource(cluster_mod.cluster)
        assert "max_k = min(max_edges + 1, len(norm_embeddings))" in source


class TestCalculateN50:
    """Test calculateN50 edge cases."""

    def test_single_contig(self):
        from COMEBin.utils import calculateN50
        assert calculateN50([1000]) == 1000

    def test_two_contigs(self):
        from COMEBin.utils import calculateN50
        assert calculateN50([1000, 500]) == 500


class TestFileLen:
    """Test file_len utility."""

    def test_empty_file(self):
        from COMEBin.utils import file_len
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False)
        tmp.close()
        try:
            result = file_len(tmp.name)
            assert result == 0
        finally:
            os.unlink(tmp.name)

    def test_non_empty_file(self):
        from COMEBin.utils import file_len
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False)
        tmp.write("line1\nline2\nline3\n")
        tmp.close()
        try:
            result = file_len(tmp.name)
            assert result == 3
        finally:
            os.unlink(tmp.name)


class TestGetLength:
    """Test get_length for FASTA files."""

    def test_empty_file(self):
        from COMEBin.utils import get_length
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        tmp.close()
        try:
            result = get_length(tmp.name)
            assert result == {}, "Empty file should return empty dict"
        finally:
            os.unlink(tmp.name)

    def test_fasta_contigs(self):
        from COMEBin.utils import get_length
        tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
        tmp.write(">contig1\nACGTACGTA\n>contig2\nACGTACGTACGT\n")
        tmp.close()
        try:
            result = get_length(tmp.name)
            assert result == {"contig1": 9, "contig2": 12}
        finally:
            os.unlink(tmp.name)
