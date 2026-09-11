# -*- coding: utf-8 -*-
#
#  One-shot script to convert chytorch rxnmap weights.pt to ONNX format.
#  Requires: torch, numpy
#  Usage: python scripts/convert_rxnmap_to_onnx.py [path/to/weights.pt] [output.onnx]
#
import sys
from importlib import import_module
from math import sqrt
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class GraphormerAttention(nn.Module):
    def __init__(self, embed_dim, num_heads, return_weights=False):
        super().__init__()
        self.embed_dim = embed_dim
        self.num_heads = num_heads
        self._scale = 1.0 / sqrt(embed_dim / num_heads)
        self.qkv_proj = nn.Linear(embed_dim, 3 * embed_dim)
        self.o_proj = nn.Linear(embed_dim, embed_dim)
        self.return_weights = return_weights

    def forward(self, x, attn_mask):
        # x: [1, N, D], attn_mask: [1, H, N, N]
        b, n, _ = x.shape
        q, k, v = self.qkv_proj(x).chunk(3, dim=-1)
        # Reshape: [1, N, H*E] -> [1, H, N, E]
        h = self.num_heads
        e = self.embed_dim // h
        q = q.view(b, n, h, e).permute(0, 2, 1, 3)  # [1, H, N, E]
        k = k.view(b, n, h, e).permute(0, 2, 3, 1)  # [1, H, E, N]
        v = v.view(b, n, h, e).permute(0, 2, 1, 3)  # [1, H, N, E]

        a = (q @ k) * self._scale + attn_mask  # [1, H, N, N]
        a = F.softmax(a, dim=-1)

        o = (a @ v).permute(0, 2, 1, 3).reshape(b, n, self.embed_dim)
        o = self.o_proj(o)

        if self.return_weights:
            return o, a.sum(dim=1) / h  # average over heads
        return o


class EncoderLayer(nn.Module):
    def __init__(self, d_model, nhead, dim_feedforward, return_weights=False):
        super().__init__()
        self.self_attn = GraphormerAttention(d_model, nhead, return_weights=return_weights)
        self.linear1 = nn.Linear(d_model, dim_feedforward)
        self.linear2 = nn.Linear(dim_feedforward, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.return_weights = return_weights

    def forward(self, x, attn_mask):
        # Post-norm transformer (norm_first=False)
        e = self.self_attn(x, attn_mask)
        if self.return_weights:
            e, a = e
        x = self.norm1(x + e)
        x = self.norm2(x + self.linear2(F.gelu(self.linear1(x))))
        if self.return_weights:
            return x, a
        return x


class RxnMapModel(nn.Module):
    """Full reaction attention model matching chytorch rxnmap architecture.

    Architecture:
    - MoleculeEncoder: atoms_encoder(121,1024) + neighbors_encoder(17,1024) + spatial_encoder(13,16)
      + 1 shared EncoderLayer applied 8 times (16 heads)
    - ReactionEncoder: role_encoder(4,1024)
      + 1 shared EncoderLayer applied 8 times (4 heads)
      Last iteration returns attention weights averaged across heads.
    """

    def __init__(self):
        super().__init__()
        d_model = 1024
        # Molecule encoder embeddings
        self.atoms_encoder = nn.Embedding(121, d_model, padding_idx=0)
        self.neighbors_encoder = nn.Embedding(17, d_model, padding_idx=0)
        self.spatial_encoder = nn.Embedding(13, 16)  # neg_inf at idx 0 handled by weights

        # Molecule encoder transformer (shared weights, applied 8 times)
        self.mol_layer = EncoderLayer(d_model, 16, 3072)

        # Role encoder
        self.role_encoder = nn.Embedding(4, d_model, padding_idx=0)

        # Reaction encoder transformer (shared weights, applied 8 times)
        # Last iteration returns attention weights
        self.rxn_layer = EncoderLayer(d_model, 4, 3072)
        self.rxn_layer_last = EncoderLayer(d_model, 4, 3072, return_weights=True)

    def forward(self, atoms, neighbors, distances, roles):
        """
        atoms: [1, N] int64
        neighbors: [1, N] int64
        distances: [1, N, N] int64
        roles: [1, N] int64

        Returns: attention matrix [N, N]
        """
        n = atoms.size(1)

        # --- Molecule encoder ---
        x = self.atoms_encoder(atoms) + self.neighbors_encoder(neighbors)  # [1, N, 1024]

        # Spatial bias (distance encoding) - same for all 8 layers (shared_attention_bias=True)
        d_mask_mol = self.spatial_encoder(distances).permute(0, 3, 1, 2)  # [1, 16, N, N]

        # Apply shared layer 8 times
        for _ in range(8):
            x = self.mol_layer(x, d_mask_mol)

        # --- Role masking and encoding ---
        # Zero out cls tokens (roles <= 1 means rxn_cls or mol_cls)
        mask = (roles > 1).unsqueeze(-1).float()  # [1, N, 1]
        x = x * mask
        x = x + self.role_encoder(roles)

        # --- Reaction encoder ---
        # Create reaction-level attention mask: -inf for padding (roles == 0)
        d_mask_rxn = torch.zeros_like(roles, dtype=torch.float32)
        d_mask_rxn = d_mask_rxn.masked_fill(roles == 0, float('-inf'))
        d_mask_rxn = d_mask_rxn.view(1, 1, 1, n).expand(1, 4, n, n)  # [1, 4, N, N]

        # Apply shared layer 7 times
        for _ in range(7):
            x = self.rxn_layer(x, d_mask_rxn)

        # Last iteration: get attention weights
        _, a = self.rxn_layer_last(x, d_mask_rxn)
        return a[0]  # [N, N]


def load_old_weights(model, weights_path):
    """Load old-format weights.pt into new model structure."""
    old = torch.load(weights_path, map_location='cpu')

    state = {}
    # Molecule encoder embeddings (rename centrality -> neighbors, spatial -> spatial)
    state['atoms_encoder.weight'] = old['molecule_encoder.atoms_encoder.weight']
    state['neighbors_encoder.weight'] = old['molecule_encoder.centrality_encoder.weight']
    state['spatial_encoder.weight'] = old['molecule_encoder.spatial_encoder.weight']

    # Molecule encoder transformer layer (rename in_proj -> qkv_proj, out_proj -> o_proj)
    state['mol_layer.self_attn.qkv_proj.weight'] = old['molecule_encoder.layer.self_attn.in_proj_weight']
    state['mol_layer.self_attn.qkv_proj.bias'] = old['molecule_encoder.layer.self_attn.in_proj_bias']
    state['mol_layer.self_attn.o_proj.weight'] = old['molecule_encoder.layer.self_attn.out_proj.weight']
    state['mol_layer.self_attn.o_proj.bias'] = old['molecule_encoder.layer.self_attn.out_proj.bias']
    state['mol_layer.linear1.weight'] = old['molecule_encoder.layer.linear1.weight']
    state['mol_layer.linear1.bias'] = old['molecule_encoder.layer.linear1.bias']
    state['mol_layer.linear2.weight'] = old['molecule_encoder.layer.linear2.weight']
    state['mol_layer.linear2.bias'] = old['molecule_encoder.layer.linear2.bias']
    state['mol_layer.norm1.weight'] = old['molecule_encoder.layer.norm1.weight']
    state['mol_layer.norm1.bias'] = old['molecule_encoder.layer.norm1.bias']
    state['mol_layer.norm2.weight'] = old['molecule_encoder.layer.norm2.weight']
    state['mol_layer.norm2.bias'] = old['molecule_encoder.layer.norm2.bias']

    # Role encoder
    state['role_encoder.weight'] = old['role_encoder.weight']

    # Reaction encoder transformer layer (both rxn_layer and rxn_layer_last share same weights)
    for prefix in ('rxn_layer', 'rxn_layer_last'):
        state[f'{prefix}.self_attn.qkv_proj.weight'] = old['layer.self_attn.in_proj_weight']
        state[f'{prefix}.self_attn.qkv_proj.bias'] = old['layer.self_attn.in_proj_bias']
        state[f'{prefix}.self_attn.o_proj.weight'] = old['layer.self_attn.out_proj.weight']
        state[f'{prefix}.self_attn.o_proj.bias'] = old['layer.self_attn.out_proj.bias']
        state[f'{prefix}.linear1.weight'] = old['layer.linear1.weight']
        state[f'{prefix}.linear1.bias'] = old['layer.linear1.bias']
        state[f'{prefix}.linear2.weight'] = old['layer.linear2.weight']
        state[f'{prefix}.linear2.bias'] = old['layer.linear2.bias']
        state[f'{prefix}.norm1.weight'] = old['layer.norm1.weight']
        state[f'{prefix}.norm1.bias'] = old['layer.norm1.bias']
        state[f'{prefix}.norm2.weight'] = old['layer.norm2.weight']
        state[f'{prefix}.norm2.bias'] = old['layer.norm2.bias']

    model.load_state_dict(state)
    return model


def main():
    # Resolved from the installed package rather than written down: no path here names a machine.
    weights_path = sys.argv[1] if len(sys.argv) > 1 else \
        str(Path(import_module('chytorch.zoo.rxnmap').__file__).parent / 'weights.pt')
    output_path = sys.argv[2] if len(sys.argv) > 2 else \
        str(Path(__file__).resolve().parent.parent / 'rxnmap.onnx')

    print(f'Loading weights from: {weights_path}')
    model = RxnMapModel()
    model = load_old_weights(model, weights_path)
    model.eval()

    # Verify model works with a dummy input
    N = 15  # example sequence length
    atoms = torch.randint(0, 121, (1, N), dtype=torch.int64)
    neighbors = torch.randint(0, 17, (1, N), dtype=torch.int64)
    distances = torch.randint(0, 13, (1, N, N), dtype=torch.int64)
    roles = torch.randint(0, 4, (1, N), dtype=torch.int64)

    with torch.no_grad():
        out = model(atoms, neighbors, distances, roles)
    print(f'Test forward pass output shape: {out.shape}')
    assert out.shape == (N, N), f'Expected ({N}, {N}), got {out.shape}'

    # Export to ONNX
    print(f'Exporting ONNX to: {output_path}')
    torch.onnx.export(
        model,
        (atoms, neighbors, distances, roles),
        output_path,
        input_names=['atoms', 'neighbors', 'distances', 'roles'],
        output_names=['attention'],
        dynamic_axes={
            'atoms': {1: 'seq_len'},
            'neighbors': {1: 'seq_len'},
            'distances': {1: 'seq_len', 2: 'seq_len'},
            'roles': {1: 'seq_len'},
            'attention': {0: 'seq_len', 1: 'seq_len'},
        },
        opset_version=17,
        do_constant_folding=True,
    )

    # Verify ONNX model
    import onnxruntime as ort
    session = ort.InferenceSession(output_path, providers=['CPUExecutionProvider'])

    # Test with different sequence length to verify dynamic axes
    N2 = 20
    atoms2 = np.random.randint(0, 121, (1, N2), dtype=np.int64)
    neighbors2 = np.random.randint(0, 17, (1, N2), dtype=np.int64)
    distances2 = np.random.randint(0, 13, (1, N2, N2), dtype=np.int64)
    roles2 = np.random.randint(0, 4, (1, N2), dtype=np.int64)

    onnx_out = session.run(None, {
        'atoms': atoms2,
        'neighbors': neighbors2,
        'distances': distances2,
        'roles': roles2,
    })[0]
    print(f'ONNX output shape with N={N2}: {onnx_out.shape}')
    assert onnx_out.shape == (N2, N2)

    # Numerical comparison
    with torch.no_grad():
        torch_out = model(
            torch.from_numpy(atoms2),
            torch.from_numpy(neighbors2),
            torch.from_numpy(distances2),
            torch.from_numpy(roles2),
        ).numpy()

    diff = np.abs(torch_out - onnx_out).max()
    print(f'Max absolute difference between torch and ONNX: {diff:.2e}')
    assert diff < 1e-4, f'Numerical mismatch: {diff}'

    import os
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f'ONNX model size: {size_mb:.1f} MB')
    print('Done!')


if __name__ == '__main__':
    main()
