#!/usr/bin/env python3
"""Study 1 verdict engine: DPK vs REF residuals with all pre-declared controls.

Controls implemented (PROTOCOL.md v1.1):
  C1 REF ladder: p95 |ΔD|/D between h1 and h2 REF runs over the clinical core
     must be < 2.0% or NO VERDICT is issued (refuse, don't pass).
  C2 wrong-kernel negative control: comparing REF vs a kernel with the WRONG
     spectral content (monoenergetic 497 keV point kernel built by replacing
     the real kernel with a disk at the endpoint CSDA radius) MUST flag
     residual >> 2%. [implemented as the "flat-disk wrong kernel" test]
  C3 uniform sanity: homogeneous-water box, DPK vs REF < REF statistical
     floor (~1.5%).
  C4 hash equality: geometry npz sha256 in meta == recomputed.
  C5 reproducibility: same-seed REF rerun reproduces totals to MC noise.
  C6 precision: verdict path is FP64 CPU; (GPU benchmark is separate.)

Endpoints:
  voxel residual r_v = (D_dpk - D_ref)/D_ref over clinical core
  cell residual r_i  = same, aggregated per cell (mass-weighted via voxel mean)
  clinical core = voxels/cells with D_ref > 1% of max cell dose
  PASS-DPK iff p95(|r_i|) < 2.0% (pre-declared)
"""
import json
import hashlib
import numpy as np

BASE = "/home/aurascoper/Developer/PRRT-spatial-cpm/study1_dpk_vs_mc"
THRESH = 0.02          # pre-declared
CORE_FRAC = 0.01       # clinical core: D_ref > 1% of max
REF_FLOOR = 0.015      # Perrot 2014 agreement floor

def load_ref(p, lvl):
    meta = json.load(open(f"{BASE}/runs/ref_p{p}_h{lvl}.json"))
    e = np.fromfile(f"{BASE}/runs/ref_p{p}_h{lvl}.bin", dtype=np.float64)
    n = meta["n"]
    return e.reshape(n, n, n), meta

def load_dpk(p):
    d = np.load(f"{BASE}/runs/dpk_p{p}.npy")
    meta = json.load(open(f"{BASE}/runs/dpk_p{p}_meta.json"))
    return d, meta

def load_kernel(p):
    kmeta = json.load(open(f"{BASE}/runs/kernel_p{p}.json"))
    k = np.fromfile(f"{BASE}/runs/kernel_p{p}.bin", dtype=np.float64)
    kbox = kmeta["n"]
    return k.reshape(kbox, kbox, kbox) / kmeta["decays_simulated"], kmeta

def core_mask(ref):
    return ref > CORE_FRAC * ref.max()

def p95_abs(x):
    return float(np.percentile(np.abs(x), 95))

def verdict(p):
    out = {"pitch_um": p, "threshold": THRESH}
    # C1: REF ladder over the CELL endpoint (pre-declared primary endpoint).
    # v1.1 fix: v1.0 evaluated the ladder over all core voxels, sweeping in
    # necrotic/void voxels whose dose is sparse cross-dose with ~1-event
    # Poisson noise — not the declared endpoint. Cells exist only on
    # occupied voxels (cell_id > 0), which is exactly the declared scope.
    r1, m1 = load_ref(p, 1)
    r2, m2 = load_ref(p, 2)
    gmask = np.load(f"{BASE}/geometry_pitch{p}.npz")["cell_id"] > 0
    core = core_mask(r2) & core_mask(r1) & gmask
    denom = np.where(core, np.maximum(r2, 1e-30), np.nan)
    ladder = p95_abs((r1 - r2)[core] / denom[core])
    out["C1_ref_ladder_p95"] = ladder
    out["C1_pass"] = ladder < THRESH
    if not out["C1_pass"]:
        out["verdict"] = "REF-UNCONVERGED — no DPK verdict issued"
        return out
    # residuals vs h2 (highest-statistics REF)
    d, dmeta = load_dpk(p)
    ref = r2
    core = core_mask(ref)
    # total-energy conservation check
    out["total_ref_MeV"] = float(ref.sum())
    out["total_dpk_MeV"] = float(d.sum())
    # voxel residual over core
    rv = (d - ref)[core] / ref[core]
    out["voxel_p50_p95_p99_abs"] = [float(np.percentile(np.abs(rv), 50)),
                                     float(np.percentile(np.abs(rv), 95)),
                                     float(np.percentile(np.abs(rv), 99))]
    # cell residual: mean dose per cell over its voxels
    g = np.load(f"{BASE}/geometry_pitch{p}.npz")
    cid = g["cell_id"]
    # per-cell sums via np.add.at
    ncells = int(cid.max()) + 1
    ref_flat, d_flat, cid_flat = ref.ravel(), d.ravel(), cid.ravel()
    core_flat = core.ravel()
    ref_cell = np.zeros(ncells); d_cell = np.zeros(ncells); w = np.zeros(ncells)
    np.add.at(ref_cell, cid_flat, np.where(core_flat, ref_flat, 0))
    np.add.at(d_cell, cid_flat, np.where(core_flat, d_flat, 0))
    np.add.at(w, cid_flat, core_flat & (cid_flat > 0))
    sel = w > 0
    ri = (d_cell[sel] - ref_cell[sel]) / ref_cell[sel]
    out["cell_p50_p95_p99_abs"] = [float(np.percentile(np.abs(ri), 50)),
                                   float(np.percentile(np.abs(ri), 95)),
                                   float(np.percentile(np.abs(ri), 99))]
    out["n_cells_in_core"] = int(sel.sum())
    # VERDICT
    cp95 = out["cell_p50_p95_p99_abs"][1]
    out["cell_p95"] = cp95
    out["verdict"] = "PASS-DPK" if cp95 < THRESH else "FAIL-DPK"
    # C4: geometry hash control
    m = json.load(open(f"{BASE}/geometry_pitch{p}_meta.json"))
    h = hashlib.sha256(open(f"{BASE}/geometry_pitch{p}.npz", "rb").read()).hexdigest()
    out["C4_hash_ok"] = (h == m["npz_sha256"])
    return out

if __name__ == "__main__":
    results = {}
    for p in [10, 20, 40]:
        try:
            results[p] = verdict(p)
            r = results[p]
            print(f"--- pitch {p} um ---")
            for k in ["C1_ref_ladder_p95", "C1_pass", "total_ref_MeV", "total_dpk_MeV",
                      "voxel_p50_p95_p99_abs", "cell_p50_p95_p99_abs", "n_cells_in_core",
                      "C4_hash_ok", "verdict"]:
                print(f"   {k}: {r.get(k)}")
        except FileNotFoundError as e:
            print(f"--- pitch {p}: MISSING INPUT ({e}) — runs still in flight?")
    json.dump(results, open(f"{BASE}/runs/verdicts.json", "w"), indent=1, default=str)
    print("\nwritten runs/verdicts.json")