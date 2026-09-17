#!/usr/bin/env pvpython
"""4D histology export: 3D space x 4 PRRT cycles, from the recorded run.

Replays the deterministic population dynamics of the representative closed
replicate (dyn seed 20261001) against the RECORDED per-cycle dose fields:
  - cycle k dose field: runs/closed_c{k}_d20261001_dose.bin (transported)
  - SF_i per occupied cell from the locked Step-2 LQ
  - deaths: Bernoulli(1-SF) with the declared dynamics seed
  - refill: well-mixed Poisson-division to capacity, drift sigma 0.10

Snapshots per cycle: occupancy (cell_id), expression e_i (voxel map),
ln e, uptake map (K1), dose (Gy), SF. Written as one .vti per cycle
(same grid, time-varying arrays) + .pvd series => scrubable 4D in ParaView.
"""
import json
import os
import numpy as np
import vtk
from vtk.util import numpy_support as vns

S1 = "/home/aurascoper/Developer/PRRT-spatial-cpm/study1_dpk_vs_mc"
S3 = "/home/aurascoper/Developer/PRRT-spatial-cpm/study3_closed_loop"
OUT = f"{S3}/paraview"
os.makedirs(OUT, exist_ok=True)

gmeta = json.load(open(f"{S1}/geometry_pitch40_meta.json"))
gd = np.load(f"{S1}/geometry_pitch40.npz")
cell_id0 = gd["cell_id"].astype(np.int32)
tissue = gd["tissue_type"].astype(np.int8)
E_T0_vox = gd["activity"].astype(np.float64)      # T0 trait map (voxel space)
n = int(gmeta["n_voxels_per_axis"])
NV = n ** 3
PITCH_UM = float(gmeta["pitch_um"])
VOX = np.flatnonzero(cell_id0 > 0).astype(np.int64)
CAP = VOX.size

verdict = json.load(open(f"{S3}/verdict.json"))
kphys = verdict["kphys_Gy_per_MeV_per_decay"]
G96 = float(verdict["declared"]["G96"])
ALPHA = float(verdict["declared"]["alpha"])
BETA = float(verdict["declared"]["beta"])
SIGMA_DIV = 0.10
T_D = 58.4
LAMBDA_DIV = np.log(2.0) / T_D
DT_DAY = 0.25
P_STEP = 1.0 - np.exp(-LAMBDA_DIV * 24.0 * DT_DAY)
N_STEPS = 224
SEED_DYN = 20261001               # representative closed replicate
SEED_DECLARED_NOTE = "replay uses recorded dose fields + declared dynamics seed"


def sf_of_D(D):
    return np.exp(-(ALPHA * D + BETA * G96 * D ** 2))


def dose_Gy_for_cycle(c):
    meta = json.load(open(f"{S3}/runs/closed_c{c}_d20261001_dose.json"))
    E = np.fromfile(f"{S3}/runs/closed_c{c}_d20261001_dose.bin", dtype=np.float64)
    E = E / meta["decays_simulated"]
    return kphys * E.reshape((n, n, n), order="C")


# ---- state: per-slot occupancy + expression; slots map to fixed voxels ----
rng = np.random.default_rng(20261001)
occ = np.ones(CAP, bool)
e = E_T0_vox.ravel(order="C")[VOX].copy()   # per-slot expression (T0)

frames = []
for cyc in range(1, 5):
    D_field = dose_Gy_for_cycle(cyc)
    D_v = D_field.ravel(order="C")[VOX[occ]]
    S_v = sf_of_D(D_v)
    # uptake map for THIS cycle: A ∝ e over occupied, K1-normalized
    uptake_field = np.zeros(NV)
    uptake_field[VOX[occ]] = e[occ]
    uptake_field *= 1.0 / uptake_field.sum()

    # snapshot BEFORE dynamics: the state the cycle's transport saw
    occ_field = np.zeros(NV, np.int32)
    occ_field[VOX[occ]] = 1
    e_field = np.zeros(NV)
    e_field[VOX[occ]] = e[occ]
    ln_e_field = np.full(NV, np.nan)
    ln_e_field[VOX[occ]] = np.log(e[occ])
    S_field = np.full(NV, np.nan)
    S_field[VOX[occ]] = 0.0
    S_field[VOX[occ]] = 0.0
    # scatter SF into voxel space for occupied cells
    tmp = np.zeros(NV)
    tmp[VOX[occ]] = S_v
    S_field = tmp
    frames.append({
        "cycle": cyc,
        "n": int(occ.sum()),
        "occ": occ_field.reshape(n, n, n),
        "e": e_field.reshape(n, n, n),
        "ln_e": ln_e_field.reshape(n, n, n),
        "uptake": (uptake_field * CAP).reshape(n, n, n),  # scaled for display
        "dose": D_field,
        "SF": S_field.reshape(n, n, n),
        "mean_ln_e": float(np.log(e[occ]).mean()),
        "mean_dose": float(D_v.mean()),
        "mean_SF": float(S_v.mean()),
    })
    print(f"cycle {cyc}: n={int(occ.sum())} mean_ln_e={np.log(e[occ]).mean():+.4f} "
          f"mean_dose={D_v.mean():.3f} Gy mean_SF={S_v.mean():.4f}")

    # ---- dynamics to next cycle (first-division death + refill)
    if cyc < 4:
        sf_slot = np.full(CAP, np.nan)
        sf_slot[occ] = S_v
        live = np.flatnonzero(occ)
        u = rng.random(live.size)
        died = live[u >= sf_slot[live]]
        occ[died] = False
        e[died] = np.nan
        deaths = int(died.size)
        divisions = 0
        for _ in range(224):
            free_idx = np.flatnonzero(~occ)
            if free_idx.size == 0:
                break
            occ_idx = np.flatnonzero(occ)
            if occ_idx.size == 0:
                break
            attempt = occ_idx[rng.random(occ_idx.size) < P_STEP]
            if attempt.size == 0:
                continue
            rng.shuffle(attempt)
            k = min(attempt.size, free_idx.size)
            tgt = free_idx[rng.permutation(free_idx.size)[:k]]
            occ[tgt] = True
            e[tgt] = e[attempt[:k]] * np.exp(SIGMA_DIV * rng.standard_normal(k))
            divisions += k
        print(f"  -> deaths {deaths}, divisions {divisions}, end n {int(occ.sum())}")

# ---- write the .vti series (one per cycle, same grid => 4D scrub in ParaView)
written = []
for fr in frames:
    img = vtk.vtkImageData()
    img.SetDimensions(n, n, n)
    mm = PITCH_UM / 1000.0
    img.SetSpacing(mm, mm, mm)
    img.SetOrigin(0.0, 0.0, 0.0)
    for name, arr, typ in [
        ("cell_id_T0", cell_id0, np.int32),
        ("tissue_type", tissue, np.int8),
        ("occupancy", fr["occ"].astype(np.uint8), np.uint8),
        ("expression_e", fr["e"], np.float64),
        ("ln_expression", fr["ln_e"], np.float64),
        ("uptake_relative", fr["uptake"], np.float64),
        ("dose_Gy", fr["dose"], np.float64),
        ("SF", fr["SF"], np.float64),
    ]:
        v = vns.numpy_to_vtk(np.ascontiguousarray(arr.ravel(order="C")), deep=1)
        v.SetName(name)
        img.GetPointData().AddArray(v)
    for nm, s in [
        ("orientation_probe_ijk", "x-fastest points, X=i, Y=j, Z=k"),
        ("units", "declared: 40 um per site; source: geometry_pitch40_meta.json"),
        ("cycle", str(fr["cycle"])),
        ("mean_ln_e", f"{fr['mean_ln_e']:+.5f}"),
        ("mean_dose_Gy", f"{fr['mean_dose']:.4f}"),
        ("provenance", "replay of dyn seed 20261001 against recorded per-cycle "
                       "transport fields (closed arm, commit 23c59ed)"),
    ]:
        sa = vtk.vtkStringArray()
        sa.SetName(nm)
        sa.InsertNextValue(s)
        img.GetFieldData().AddArray(sa)
    stem = f"{OUT}/hist4d_c{fr['cycle']}.vti"
    w = vtk.vtkXMLImageDataWriter()
    w.SetFileName(stem)
    w.SetInputData(img)
    if not w.Write():
        raise RuntimeError(f"write failed {stem}")
    written.append((fr["cycle"], stem))
    print("wrote", stem)

pvd = f"{OUT}/histology_4d.pvd"
with open(pvd, "w") as f:
    f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n')
    for cyc, stem in written:
        f.write(f'    <DataSet part="{cyc}" file="{os.path.basename(stem)}"/>\n')
    f.write("  </Collection>\n</VTKFile>\n")
print("wrote", pvd)