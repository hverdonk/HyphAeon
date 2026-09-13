"""
aeon-core/tests/test_splits.py
-----------------------------
Unit tests for aeon_core.splits: compute_fused_affinity_matrix.
"""

import numpy as np
import pytest

from aeon_core.splits import compute_fused_affinity_matrix


class TestFusedAffinityMatrix:
    """Unit tests for the multiplicative affinity tensor fusion."""

    def test_symmetry_and_diagonal_zeroing(self):
        n = 8
        rng = np.random.RandomState(42)
        attn = rng.uniform(0.01, 0.5, size=(n, n))
        mds = rng.normal(0, 1, size=(n, 4))
        taxa_repr = rng.normal(0, 1, size=(n, 16))

        A = compute_fused_affinity_matrix(attn, mds, taxa_repr)

        assert A.shape == (n, n)
        assert np.allclose(A, A.T, atol=1e-6)
        assert np.all(np.diag(A) == 0.0)
        assert np.all(A >= 0.0)
        assert np.all(np.isfinite(A))

    def test_without_taxa_repr(self):
        n = 6
        rng = np.random.RandomState(42)
        attn = rng.uniform(0.01, 0.5, size=(n, n))
        mds = rng.normal(0, 1, size=(n, 4))

        A = compute_fused_affinity_matrix(attn, mds, taxa_repr=None)
        assert A.shape == (n, n)
        assert np.allclose(A, A.T, atol=1e-6)
        assert np.all(np.diag(A) == 0.0)
        assert np.all(A >= 0.0)

    def test_degenerate_zero_distance_mds(self):
        n = 5
        attn = np.ones((n, n), dtype=np.float64)
        mds = np.zeros((n, 4), dtype=np.float64)

        A = compute_fused_affinity_matrix(attn, mds)
        assert A.shape == (n, n)
        assert np.all(np.isfinite(A))
        assert np.all(np.diag(A) == 0.0)
