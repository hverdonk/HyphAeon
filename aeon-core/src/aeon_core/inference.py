"""
aeon_core/inference.py
---------------------
Shared inference helpers: device selection, model loading, alignment
preparation, and adaptive batch sizing.

These are used by both hyphaeon and chronaeon packages.
"""

import os
from typing import Optional
import numpy as np
import torch

from .model import PhyloAxialTransformer
from .weights import resolve_weights_path, load_arch_config, load_weights
from .dataset import load_alignment_and_tree


def get_device_memory_budget(device: torch.device) -> float:
    """
    Dynamically probes accelerator VRAM, Unified Memory, or host RAM capacity
    across CUDA, MPS, TPU (XLA), and CPU to determine an optimal per-layer
    attention allocation budget (in bytes).
    """
    dev_type = device.type
    
    if dev_type == 'cuda' and torch.cuda.is_available():
        try:
            total_vram = torch.cuda.get_device_properties(device).total_memory
            # Allocate up to 20% of VRAM for peak single-layer attention tensor (bounded between 1.5 GB and 8.0 GB)
            return float(min(8.0e9, max(1.5e9, total_vram * 0.20)))
        except Exception:
            return 3.0e9
            
    elif dev_type == 'mps':
        try:
            # Query macOS unified RAM capacity
            sys_ram = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
            # 10% of unified memory (bounded between 1.5 GB and 6.0 GB)
            return float(min(6.0e9, max(1.5e9, sys_ram * 0.10)))
        except Exception:
            return 2.0e9
            
    elif dev_type == 'xla':
        # TPU v2/v3/v4/v5e typically feature 16GB-32GB HBM per core
        return 4.0e9
        
    else: # CPU
        try:
            sys_ram = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
            return float(min(4.0e9, max(1.0e9, sys_ram * 0.10)))
        except Exception:
            return 1.5e9


def compute_adaptive_safe_batch_size(
    n_taxa: int, 
    user_batch_size: Optional[int] = None,
    device: Optional[torch.device] = None
) -> int:
    """
    Computes an optimal, memory-safe batch size dynamically adapted to:
      1. Sequence depth N (quadratic scaling: 48 * N^2 bytes per sequence)
      2. Hardware accelerator capacity (CUDA, MPS, TPU/XLA, or CPU)
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else ('mps' if torch.backends.mps.is_available() else 'cpu'))
        
    budget_bytes = get_device_memory_budget(device)
    # Memory per sequence = 12 heads * N^2 * 4 bytes = 48 * N^2 bytes
    safe_max_b = max(1, int(budget_bytes / (48.0 * max(1, n_taxa) ** 2)))
    
    # Structure safe batches in clean multiples of 19 (for 19 AA DMS passes)
    if n_taxa >= 600:
        base_cap = 38 if budget_bytes <= 3.0e9 else (76 if budget_bytes <= 6.0e9 else 152)
        safe_b = min(base_cap, safe_max_b)
    elif n_taxa >= 400:
        base_cap = 57 if budget_bytes <= 3.0e9 else (114 if budget_bytes <= 6.0e9 else 228)
        safe_b = min(base_cap, safe_max_b)
    elif n_taxa >= 250:
        base_cap = 76 if budget_bytes <= 3.0e9 else (152 if budget_bytes <= 6.0e9 else 256)
        safe_b = min(base_cap, safe_max_b)
    elif n_taxa >= 100:
        safe_b = min(128 if budget_bytes <= 3.0e9 else 256, safe_max_b)
    else:
        safe_b = min(256 if budget_bytes <= 3.0e9 else 512, safe_max_b)
        
    safe_b = max(1, safe_b)
    
    if user_batch_size is not None and user_batch_size > 0:
        if user_batch_size > safe_b * 2:
            print(f"[!] Warning: Requested batch size {user_batch_size} exceeds hardware safety threshold ({safe_b}) for N={n_taxa} on {device.type.upper()}.", flush=True)
            print(f"    Automatically capping batch size to {safe_b} to guarantee stable execution.", flush=True)
            return safe_b
        return user_batch_size
        
    return safe_b


def get_device(cpu: bool = False) -> torch.device:
    """Select the best available hardware device.

    CUDA > MPS > CPU, unless cpu=True forces CPU.
    """
    if cpu:
        return torch.device('cpu')
    if torch.cuda.is_available():
        return torch.device('cuda')
    if torch.backends.mps.is_available():
        return torch.device('mps')
    return torch.device('cpu')


def load_model(weights=None, variant=None, device=None, strict=False):
    """Load a PhyloAxialTransformer from weights path or HuggingFace variant.

    Returns an eval-mode model on the specified device.
    """
    if device is None:
        device = get_device()

    weights_path = resolve_weights_path(weights=weights, variant=variant)
    config = load_arch_config(weights=weights, variant=variant)
    model = PhyloAxialTransformer(
        embed_dim=config['embed_dim'],
        num_layers=config['num_layers'],
        num_heads=config['num_heads'],
        window_size=config['window_size'],
    ).to(device)
    state_dict = load_weights(weights=weights_path, variant=variant, map_location=device)
    model.load_state_dict(state_dict, strict=strict)
    model.eval()
    return model


def prepare_alignment(alignment_path, tree_path=None, model=None, device=None,
                      max_species=None, prune_duplicates=True, use_tn93=False):
    """Load an alignment + tree and precompute the tree attention cache.

    If use_tn93=True, skips the phylogenetic tree and computes pairwise
    distances directly from the alignment using TN93.

    Returns (c, a, d, z, inv, taxa, L, tree_cache).
    If model is None, tree_cache will be None.
    """
    c, a, d, z, inv, taxa, L = load_alignment_and_tree(
        alignment_path, tree_path, max_species=max_species,
        prune_duplicates=prune_duplicates, use_tn93=use_tn93
    )
    tree_cache = None
    if model is not None:
        tree_cache = model.precompute_tree_cache(d.to(device), z.to(device))
    return c, a, d, z, inv, taxa, L, tree_cache
