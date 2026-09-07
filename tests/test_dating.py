"""Unit tests for HyphAeon dating module and non-linear clock models."""
import numpy as np
import pytest

from hyphaeon.dating import (
    run_ols_dating,
    run_pgls_dating,
    run_powerlaw_clock_dating,
    compute_rcs_basis,
    run_restricted_spline_clock_dating,
    run_mrca_dating,
)


class TestClockDatingModels:
    def test_ols_clock_linear(self):
        np.random.seed(42)
        true_t0 = 1950.0
        true_mu = 0.002
        times = np.linspace(1970, 2020, 50)
        noise = np.random.normal(0, 0.002, size=len(times))
        dists = true_mu * (times - true_t0) + noise

        res = run_ols_dating(times, dists, n_boot=100)
        assert abs(res['t_mrca'] - true_t0) < 5.0
        assert abs(res['mu'] - true_mu) < 0.0005
        assert res['r2'] > 0.90

    def test_rcs_basis_ancestral_linearity(self):
        """Verify Harrell RCS basis vanishes identically for t <= t_min, ensuring strictly linear ancestral extrapolation."""
        knots = np.array([1980.0, 2000.0, 2020.0])
        ancestral_times = np.array([1900.0, 1950.0, 1979.9, 1980.0])
        B, dB = compute_rcs_basis(ancestral_times, knots)
        assert np.all(B == 0.0)
        assert np.all(dB == 0.0)

    def test_restricted_spline_on_linear_synthetic(self):
        """When data is strictly linear, restricted spline should fail to reject linearity and delta_AIC < 0."""
        np.random.seed(42)
        true_t0 = 1960.0
        true_mu = 0.0015
        times = np.linspace(1970, 2020, 60)
        noise = np.random.normal(0, 0.001, size=len(times))
        dists = true_mu * (times - true_t0) + noise

        res = run_restricted_spline_clock_dating(times, dists, n_boot=50)
        assert abs(res['t_mrca'] - true_t0) < 5.0
        assert not res['is_nonlinear_preferred']

    def test_restricted_spline_on_decelerating_synthetic(self):
        """When data has genuine rate deceleration across the observation window,

        restricted spline should detect it (p < 0.05, delta_AIC > 2) without boundary collapse.
        """
        np.random.seed(123)
        times = np.linspace(1980, 2020, 80)
        knots = np.array([1980.0, 2000.0, 2015.0])
        B, _ = compute_rcs_basis(times, knots)
        # beta_0 = -1.98, beta_1 = 0.0010 (ancestral rate), beta_2 = -0.0005 (deceleration)
        true_dists = -1.98 + 0.0010 * times - 0.0005 * B[:, 0]
        noise = np.random.normal(0, 0.0002, size=len(times))
        dists = true_dists + noise

        res = run_restricted_spline_clock_dating(times, dists, n_boot=50)
        assert res['is_nonlinear_preferred']
        assert res['p_f_test'] < 0.01
        assert res['delta_aic'] > 2.0
        assert res['rate_recent'] < res['rate_ancestral']
        # Ancestral extrapolation is linear and stable
        assert 1970.0 < res['t_mrca'] < 1982.0

    def test_run_mrca_dating_auto_selection(self, tmp_path):
        """Test run_mrca_dating end-to-end with mock alignment and auto clock selection."""
        fasta_path = tmp_path / 'mock_dated.fasta'
        fasta_content = (
            ">seq1_1980\nATGGCC\n"
            ">seq2_1990\nATGGCA\n"
            ">seq3_2000\nATGGTA\n"
            ">seq4_2010\nTTGGTA\n"
            ">seq5_2020\nTTGGTT\n"
        )
        fasta_path.write_text(fasta_content)

        res = run_mrca_dating(
            alignment_path=str(fasta_path),
            use_tn93=True,
            clock_model='auto',
            method='ols',
            n_bootstrap=50,
        )

        assert 'ols' in res
        assert 'spline' in res
        assert res['clock_model'] == 'auto'
        assert 'selected_clock' in res
        assert len(res['taxa_records']) == 5

    def test_time_decay_consensus(self):
        """Test that time-decay consensus properly upweights early ancestral isolates."""
        from hyphaeon.dating import generate_time_decay_consensus_sequence
        seq_dict = {
            'seq1_1920': 'AAAAAA',
            'seq2_2020': 'TTTTTT',
            'seq3_2021': 'TTTTTT',
            'seq4_2022': 'TTTTTT',
        }
        dates_map = {'seq1_1920': 1920.0, 'seq2_2020': 2020.0, 'seq3_2021': 2021.0, 'seq4_2022': 2022.0}
        # In unweighted consensus, 'T' is 3/4 (75%), so unweighted would produce 'TTTTTT'.
        # But with time-decay (gamma=0.05, delta_t=100 -> weight of 1920 is exp(0)=1.0 vs ~exp(-5)=0.0067),
        # 1920 dominates and produces 'AAAAAA'.
        con_seq, eff_g = generate_time_decay_consensus_sequence(seq_dict, dates_map, gamma=0.05)
        assert con_seq == 'AAAAAA'
        assert eff_g == 0.05
