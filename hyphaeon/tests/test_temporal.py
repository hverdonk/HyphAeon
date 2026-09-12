"""Unit tests for temporal selection dynamics in hyphaeon.temporal."""
import numpy as np
import pandas as pd
import pytest
import torch

from hyphaeon.temporal import (
    parse_temporal_metadata,
    infer_root_sequence,
)


class TestTemporalSelectionDynamics:
    def test_parse_temporal_metadata_table(self, tmp_path):
        csv_file = tmp_path / "meta.csv"
        df = pd.DataFrame({
            "strain": ["seq_A", "seq_B", "seq_C"],
            "collection_date": ["2020-03-15", "2021-06-20", "2022.5"]
        })
        df.to_csv(csv_file, index=False)

        dates = parse_temporal_metadata(csv_file)
        assert len(dates) == 3
        assert abs(dates["seq_A"] - 2020.20) < 0.02
        assert abs(dates["seq_B"] - 2021.46) < 0.02
        assert dates["seq_C"] == 2022.5

    def test_infer_root_sequence(self):
        a = torch.tensor([
            [[0], [0], [0], [5]],
            [[1], [1], [3], [3]],
            [[2], [2], [2], [2]],
        ], dtype=torch.long)
        taxa = ["t1", "t2", "t3", "t4"]
        taxa_dates = np.array([2020.0, 2020.1, 2022.0, 2022.5])

        root_aas, root_idx = infer_root_sequence(a, taxa, taxa_dates=taxa_dates)
        assert len(root_aas) == 3
        assert root_idx[0] == 0
        assert root_idx[1] == 1
        assert root_idx[2] == 2
