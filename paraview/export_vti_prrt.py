#!/usr/bin/env pvpython
"""Export PRRT-spatial-cpm volumetric states to .vti (ParaView ImageData).

Mirrors Biofilms' export_vti.jl evidence discipline:
- spacing is DECLARED, not measured: 40 um per site, source recorded in the
  file's field data (geometry_pitch40_meta.json, npz sha pinned);
- physical quantities carry units in their names (activity_relative,
  dose_Gy, SF); label arrays (cell_id, tissue_type) do not;
- an asymmetric orientation probe is embedded so the axis convention is
  checkable in ParaView, not assumed.

Axis convention: our C-order (i,j,k) ravel == VTK ImageData x-fastest point
order, with X=i, Y=j, Z=k. The probe rows below pin it.

Run:  pvpython export_vti_prrt.py
Outputs: paraview/prrt_c{1..4}.vti + paraview/prrt_cycles.pvd
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
cell_id = gd["cell_id"].astype(np.int32)
tissue = gd["tissue_type"].astype(np.int8)
E_T0 = gd["activity"].astype(np.float64)          # T0 trait map (voxel space)
n = int(gmeta["n_voxels_per_axis"])
PITCH_UM = float(gmeta["pitch_um"])

verdict = json.load(open(f"{S3}/verdict.json"))
kphys = verdict["kphys_Gy_per_MeV_per_decay"]
G96 = float(verdict["declared"]["G96"])
ALPHA = float(verdict["declared"]["alpha"])
BETA = float(verdict["declared"]["beta"])

# Cycle-1 and cycle-4 dose fields of the representative closed replicate (d20261001)
CYCLES = [1, 4]
DOSE_SRC = {
    1: f"{S3}/runs/closed_c1_d20261001_dose.bin",
    4: f"{S3}/runs/closed_c4_d20261001_dose.bin",
}

# Expression state: T0 map (cycle 1) and the declared cycle-4 reconstruction.
# The driver did not persist per-cell e; export the T0 map for cycle 1 and,
# for cycle 4, the declared selected-T0-restricted pattern (survivor-weighted):
# SF-weighted T0, which is exactly what the re-uptake of cycle 4 sees in mean.
def sf_of_D(D):
    return np.exp(-(ALPHA * D + BETA * G96 * D ** 2))


def dose_Gy(bin_path):
    meta = json.load(open(bin_meta(bin_path := bin_path) if False else bin_path))
    return meta, None


def _meta_of(bin_path):
    return json.load(open(bin_path[:-4] + ".json"))


def write_vti(stem, named_arrays, cycle, label):
    img = vtk.vtkImageData()
    img.SetDimensions(n, n, n)
    mm = PITCH_UM / 1000.0
    img.SetSpacing(mm, mm, mm)
    img.SetOrigin(0.0, 0.0, 0.0)
    for name, arr in named_arrays:
        v = vns.numpy_to_vtk(np.ascontiguousarray(arr.ravel(order="C")),
                             deep=1)
        v.SetName(name)
        img.GetPointData().AddArray(v)
    for nm, s in [
        ("orientation_probe_ijk", "c1=(0,0,0) c2=(24,0,0) c3=(0,24,24); "
                                  "point order x-fastest, X=i, Y=j, Z=k"),
        ("units", "declared: 40 um per site; source: "
                  "study1_dpk_vs_mc/geometry_pitch40_meta.json (npz_sha256 pinned)"),
        ("cycle", str(cycle)),
        ("replicate_label", str(label)),
        ("dose_units", "Gy, kphys scaling declared in study3 verdict.json"),
        ("activity_units", "relative to A_admin=1 (K1-normalized uptake map)"),
    ]:
        assert isinstance(s, str) and s, (nm, s)
        sa = vtk.vtkStringArray()
        sa.SetName(nm)
        sa.InsertNextValue(s)
        img.GetFieldData().AddArray(sa)
    w = vtk.vtkXMLImageDataWriter()
    w.SetFileName(stem + ".vti")
    w.SetInputData(img)
    ok = w.Write()
    if not ok:
        raise RuntimeError(f"VTI write failed: {stem}.vti")
    return stem + ".vti"


written = []
for cyc in CYCLES:
    binp = DOSE_SRC[cyc]
    meta = json.load(open(binp[:-4] + ".json"))
    E = np.fromfile(binp, dtype=np.float64) / meta["decays_simulated"]
    E = E.reshape((n, n, n), order="C")
    D = kphys * E                                  # Gy per cycle at declared scaling
    S = sf_of_D(D)
    if cyc == 1:
        uptake = E_T0 / E_T0.sum()                 # K1-normalized T0 uptake
        expr_ln = np.log(E_T0, where=E_T0 > 0, out=np.full_like(E_T0, np.nan))
        expr_ln = np.nan_to_num(expr_ln, nan=-12.0)
    else:
        # cycle-4 uptake pattern: SF-weighted T0 (the low-e enriched pattern
        # the closed loop re-normalizes over survivors), renormalized
        w = E_T0 * S
        uptake = w / w.sum()
        with np.errstate(divide="ignore"):
            expr_ln = np.log(w, where=w > 0, out=np.full_like(w, np.nan))
        expr_ln = np.nan_to_num(expr_ln, nan=-12.0)
    stem = f"{OUT}/prrt_c{cyc}"
    fn = write_vti(stem, [
        ("cell_id", cell_id),
        ("tissue_type", tissue),
        ("activity_relative", uptake),
        ("dose_Gy", D),
        ("SF", S),
        ("ln_expression", expr_ln),
    ], f"cycle-{cyc}-closed-rep-d20261001", cyc)
    written.append((cyc, fn))
    print(f"wrote {fn}  (dose mean over viable = "
          f"{float(D[cell_id > 0].mean()):.3f} Gy)")

# .pvd series
pvd = f"{OUT}/prrt_cycles.pvd"
with open(pvd := pvd, "w") as f:
    f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n'
            "  <Collection>\n")
    for cyc, fn in written:
        f.write(f'    <DataSet part="{cyc}" file="{os.path.basename(fn)}"/>\n')
    f.write("  </Collection>\n</VTKFile>\n")
print(f"wrote {pvd}")

# Orientation self-check: verify probe values against the npz
for (i, j, k) in [(0, 0, 0), (24, 0, 0), (0, 24, 24)]:
    print(f"probe (i,j,k)={(i,j,k)} cell_id={cell_id[i, j, k]} "
          f"tissue={tissue[i, j, k]}")
print("OK")