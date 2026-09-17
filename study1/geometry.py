#!/usr/bin/env python3
"""Synthetic histology-representative CPM snapshot + activity map generator.

One geometry, fixed for the whole study (fixed-object ladder rule): the
snapshot is generated ONCE per pitch with declared seeds, hashed, and both
arms must consume the identical hashed files.

Structure (declared, histology-motivated volume fractions for a treated
NET lesion voxel):
  - viable tumor cell belt (outer shell, ~60% of ROI volume)
  - necrotic core (~20%)
  - stromal bands between (~20%)
Cells are placed on a slightly jittered lattice inside the viable/stroma
regions; activity attaches per-cell via the lognormal expression draw
(preprint Eq. 4 setup), scaled to a declared mean activity concentration.

Outputs per pitch: geometry.npz (cell_id field, cell types, per-cell
expression, activity map), geometry_meta.json with sha256s.
"""
import json
import hashlib
import numpy as np

PITCHS_UM = [10.0, 20.0, 40.0]
ROI_MM = 1.0
SEED_GEOM = 20260916
SEED_EXPR = 20260917

# declared cell-level parameters
CELL_DIAM_UM = 12.0            # NET cells ~10-15 um
VOLUME_FRACT = {"viable": 0.60, "necrotic": 0.20, "stroma": 0.20}
EXPR_MU = 0.0                  # lognormal(mu, sigma) in arbitrary units
EXPR_SIGMA = 0.6                # declared spread (preprint: declared prior)
MEAN_ACTIVITY_BQ_PER_VOXEL = 1.0  # arbitrary scale; both arms consume same map

def make_geometry(pitch_um):
    rng_geom = np.random.default_rng(SEED_GEOM)
    n = int(round(ROI_MM * 1000.0 / pitch_um))
    xs = (np.arange(n) + 0.5) * pitch_um / 1000.0  # mm
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt((X - 0.5)**2 + (Y - 0.5)**2 + (Z - 0.5)**2)

    # necrotic core: ball of ~ the necrotic volume fraction
    r_nec = 0.5 * (VOLUME_FRACT["necrotic"] ** (1/3))
    necrotic = r < r_nec
    # stromal band: shell around the core
    r_stroma_in, r_stroma_out = r_nec, r_nec + 0.12
    stroma = (~necrotic) & (r < r_stroma_out)
    # viable: everything else
    viable = (~necrotic) & (~stroma)
    tissue_type = np.zeros((n, n, n), dtype=np.int8)  # 0 void(unused),1 viable,2 necrotic,3 stroma
    tissue_type[viable] = 1
    tissue_type[necrotic] = 2
    tissue_type[stroma] = 3

    # cells on jittered lattice where tissue is viable or stromal
    cell_pitch_um = CELL_DIAM_UM * 1.15
    step = max(1, int(round(cell_pitch_um / pitch_um)))
    cell_id = np.zeros((n, n, n), dtype=np.int32)  # 0 = empty (declared semantic)
    next_id = 1
    for i in range(0, n, step):
        for j in range(0, n, step):
            for k in range(0, n, step):
                if tissue_type[i, j, k] in (1, 3):
                    cell_id[i, j, k] = next_id
                    next_id += 1
    n_cells = next_id - 1

    # per-cell expression + activity (lognormal draw)
    rng_expr = np.random.default_rng(SEED_EXPR)
    e = rng_expr.lognormal(EXPR_MU, EXPR_SIGMA, size=n_cells)
    # map cell expression to voxel activity
    activity = np.zeros((n, n, n))
    flat_ids = cell_id.ravel()
    flat_act = activity.ravel()
    np.add.at(flat_act, np.nonzero(flat_ids)[0], e[flat_ids[flat_ids > 0] - 1])
    # rescale to declared mean per occupied voxel
    occ = flat_ids > 0
    flat_act[occ] *= MEAN_ACTIVITY_BQ_PER_VOXEL / flat_act[occ].mean()

    meta = {
        "pitch_um": pitch_um,
        "n_voxels_per_axis": n,
        "roi_mm": ROI_MM,
        "seed_geometry": SEED_GEOM,
        "seed_expression": SEED_EXPR,
        "cell_diam_um": CELL_DIAM_UM,
        "cell_lattice_step_voxels": step,
        "volume_fractions_declared": VOLUME_FRACT,
        "volume_fractions_realized": {
            "viable": float((tissue_type == 1).mean()),
            "necrotic": float((tissue_type == 2).mean()),
            "stroma": float((tissue_type == 3).mean()),
        },
        "n_cells": int(n_cells),
        "expr_lognormal": {"mu": EXPR_MU, "sigma": EXPR_SIGMA},
        "mean_activity_bq_per_voxel": MEAN_ACTIVITY_BQ_PER_VOXEL,
        "activity_cv": float(np.std(flat_act[occ]) / np.mean(flat_act[occ])),
    }
    return cell_id, tissue_type, activity.reshape(n, n, n), meta

if __name__ == "__main__":
    for p in PITCHS_UM:
        cell_id, tissue_type, activity, meta = make_geometry(p)
        tag = f"pitch{int(p)}"
        npz = f"geometry_{tag}.npz"
        np.savez_compressed(npz, cell_id=cell_id, tissue_type=tissue_type, activity=activity)
        h = hashlib.sha256(open(npz, "rb").read()).hexdigest()
        meta["npz_sha256"] = h
        with open(f"geometry_{tag}_meta.json", "w") as f:
            json.dump(meta, f, indent=1)
        print(f"[{tag}] n={meta['n_voxels_per_axis']}^3 cells={meta['n_cells']} "
              f"CV={meta['activity_cv']:.3f} sha256={h[:16]}...")
    print("OK")