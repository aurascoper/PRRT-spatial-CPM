#!/usr/bin/env python3
"""Study 3 figure + state export for the public repo.

Produces, paths relative to the repository root:
  data/cycle4_state.npz                   — cycle-1 and cycle-4 dose fields, T0 map
  figures/activity_cycle1_vs_cycle4.png   — activity map, cycle 1 vs cycle 4
  figures/trajectory.png                  — population + drift trajectory, both arms
  figures/dose_cv_compression.png         — dose CV per cycle, both arms

Reads study3_closed_loop/runs/closed_c{1,4}_d20261001_dose.bin and
closed_c4_d20261001.bin, the representative closed-loop replicate.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S3 = os.path.dirname(os.path.abspath(__file__))      # study3_closed_loop/
DATA = os.path.join(S3, "..", "data")                # hashed geometry npz, state export
FIG = os.path.join(S3, "..", "figures")
RUNS = f"{S3}/runs"
os.makedirs(FIG, exist_ok=True)
os.makedirs(RUNS, exist_ok=True)

verdict = json.load(open(f"{S3}/verdict.json"))
gdat = np.load(f"{DATA}/geometry_pitch40.npz")
cell_id = gdat["cell_id"]
n = 25
VOX = np.flatnonzero(cell_id > 0).astype(np.int64)
CAP = VOX.size
E_T0 = gdat["activity"].astype(np.float64).ravel(order="C")

# ---- export cycle-4 state of closed rep0 (declared representative run) ------
# e_state: reconstruct from the driver is not persisted per-cycle; export the
# declared T0 map + verdict numbers instead, and save the c4 dose field.
def dose_per_decay(tag):
    m = json.load(open(f"{RUNS}/{tag}_dose.json"))
    return np.fromfile(f"{RUNS}/{tag}_dose.bin", dtype=np.float64) / m["decays_simulated"], m
d4, m4 = dose_per_decay("closed_c4_d20261001")
d1, m1 = dose_per_decay("closed_c1_d20261001")
kphys = verdict["kphys_Gy_per_MeV_per_decay"]
mean_gy = kphys * float(d1.reshape(n, n, n)[cell_id > 0].mean())
assert abs(mean_gy - 10.0) < 1e-6, f"cycle-1 viable mean dose {mean_gy} Gy, declared 10 Gy (issue #9)"
np.savez_compressed(f"{DATA}/cycle4_state.npz",
                    cell_id=cell_id,
                    dose_cycle1_per_decay=d1.reshape(n, n, n),
                    dose_cycle4_per_decay=d4.reshape(n, n, n),
                    e_T0_voxel=E_T0.reshape(n, n, n),
                    meta=json.dumps({
                        "geometry_sha256": "a6883b84a577dd5db888c9d4292cf68897a68d840b64a427061bc70a2d60a170",
                        "dose_units": f"MeV/voxel per decay; {m1['decays_simulated']} decays per cycle; "
                                      f"seeds cycle 1 {m1['seed']}, cycle 4 {m4['seed']}",
                        "kphys_Gy_per_MeV_per_decay": verdict["kphys_Gy_per_MeV_per_decay"],
                        "alpha": 0.24, "beta": 0.06, "G96": 0.044068,
                        "ln_e_drift_closed_mean": verdict["delta_ln_e_closed_mean"],
                        "ln_e_drift_open_mean": verdict["delta_ln_e_open_mean"],
                    }))
print(f"exported data/cycle4_state.npz (cycle-1 viable mean {mean_gy:.4f} Gy)")

# ---- Figure 1: expression-proportional activity, cycle 1 vs cycle 4 shape ---
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
sl = 12  # central slice
a1 = gdat["activity"][sl]          # cycle-1 pattern ∝ e_T0
a4 = np.fromfile(f"{RUNS}/closed_c4_d20261001.bin", dtype=np.float64,
                 offset=n**3 * 4).reshape(n, n, n)[sl]
vmax = max(a1.max(), a4.max())
for ax, a, ttl in [(axes[0], a1, "Cycle 1: $A_j \\propto e_i$ (T0 map)"),
                   (axes[1], a4, "Cycle 4: re-normalized uptake")]:
    im = ax.imshow(a, cmap="inferno", vmin=0, vmax=vmax, origin="lower")
    ax.set_title(ttl)
    ax.set_xlabel("y (voxel, 40 µm)")
    ax.set_ylabel("x (voxel, 40 µm)")
fig.colorbar(im, ax=axes, label="relative activity", shrink=0.85)
fig.savefig(f"{FIG}/activity_cycle1_vs_cycle4.png", dpi=160)
plt.close(fig)
print("wrote figures/activity_cycle1_vs_cycle4.png")

# ---- Figure 2: trajectory + drift -------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                               constrained_layout=True)
cyc = [1, 2, 3, 4]


def rep_mean(stats_list, key):
    """Mean over the 5 replicates of the per-cycle value."""
    return [float(np.mean([st[c][key] for st in stats_list]))
            for c in range(4)]


for arm_key, nm, col in [("closed_stats", "closed loop", "#1f77b4"),
                         ("open_stats", "open loop (frozen)", "#d62728")]:
    stl = verdict[arm_key]
    ax1.plot(cyc, rep_mean(stl, "mean_ln_e_start"), "o-", color=col,
             label=nm)
ax1.set_xlabel("cycle")
ax1.set_ylabel("population mean of ln $e_i$ at cycle start")
ax1.set_title("SSTR2 expression drift under 4-cycle PRRT")
ax1.legend()
ax1.grid(alpha=0.3)
surv_closed = [[st[c]["n_start"] - st[c]["n_deaths"] for c in range(4)]
               for st in verdict["closed_stats"]]
surv_open = [[st[c]["n_start"] - st[c]["n_deaths"] for c in range(4)]
             for st in verdict["open_stats"]]
ax2.plot(cyc, [float(np.mean([r[c] for r in surv_closed])) for c in range(4)],
         "s-", color="#1f77b4", label="closed loop (5-seed mean)")
ax2.plot(cyc, [float(np.mean([r[c] for r in surv_open])) for c in range(4)],
         "s--", color="#d62728", label="open loop (5-seed mean)")
ax2.axhline(CAP * 0.0696, color="gray", ls=":", lw=1,
            label="uniform-dose SF(mean) level")
ax2.set_xlabel("cycle")
ax2.set_ylabel("per-cycle survivors (of 13,978)")
ax2.set_title("per-cycle kill under fixed administered activity")
ax2.legend()
ax2.grid(alpha=0.3)
fig.savefig(f"{FIG}/trajectory.png", dpi=160)
print("wrote figures/trajectory.png")

# ---- Figure 3: dose-shape compression (the Jensen mechanism) -----------------
fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
for arm_key, nm, col in [("closed_stats", "closed", "#1f77b4"),
                         ("open_stats", "open", "#d62728")]:
    ax.plot(cyc, rep_mean(verdict[arm_key], "dose_cv"), "o-", color=col,
            label=nm)
ax.set_xlabel("cycle")
ax.set_ylabel("dose CV over viable cells")
ax.set_title("selection compresses the dose map (Jensen excess shrinks)")
ax.legend()
ax.grid(alpha=0.3)
fig.savefig(f"{FIG}/dose_cv_compression.png", dpi=160)
print("wrote figures/dose_cv_compression.png")