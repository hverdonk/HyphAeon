"""epistaeon/io.py -- result serialisation."""

import json
from pathlib import Path
from typing import Any, Dict, Sequence

import numpy as np


def _plain(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return [_plain(x) for x in obj.tolist()]
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(x) for x in obj]
    if isinstance(obj, float) and (obj != obj or obj in (float("inf"), float("-inf"))):
        return None            # JSON has no NaN/Infinity
    return obj


def write_json(path: str, payload: Dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(_plain(payload), indent=2))


def write_matrix_csv(path: str, matrix: np.ndarray, names: Sequence[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ["," + ",".join(names)]
    for i, row in enumerate(np.asarray(matrix)):
        cells = ["" if v != v else f"{float(v):.6g}" for v in row]
        lines.append(names[i] + "," + ",".join(cells))
    p.write_text("\n".join(lines) + "\n")
