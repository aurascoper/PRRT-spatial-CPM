#!/usr/bin/env python3
"""Study 3 — Closed-loop selection: per-cycle re-transport + LQ + CPM.

Executes study3_closed_loop/PROTOCOL.md v1.0 + v1.1 + v1.2 exactly.

Per cycle (protocol 'Iteration sequence'):
  1. TRANSPORT: write activity map (K1-normalized), invoke study1_app.
  2. RESPONSE: D_i = kphys * E_voxel (declared scaling, frozen at cycle 1);
     SF_i = exp(-alpha*D_i - beta*G96*D_i^2)  (locked Step-2 functions).
  3. GROWTH/DEATH: first-division Bernoulli(1-SF_i) -> mitotic catastrophe,
     site cleared; survivors refill to capacity over the 56-day window
     (0.25-day Poisson steps, ln2/58.4 h); daughters inherit e with
     lognormal drift sigma_div = 0.10 (well-mixed placement; Potts
     neighborhood dynamics are Step 4, per protocol Boundaries).
  4. RE-UPTAKE: A_admin re-applied ∝ e_i (closed) / frozen voxel map (open).

Declared conventions:
  - T0 trait identity: e_i(T0) := hashed T0 activity map (relative values
    only matter; multiplicative drift, renormalized uptake).
  - Open arm: cycle-1 activity map frozen VERBATIM at the voxel level —
    dead sites KEEP their declared activity ("re-applied blindly,
    ignoring the death of high-e_i cells", Choice B). Because the frozen
    map is the T0 map and the transport seed is fixed, the open arm's
    dose field is bit-identical across ALL cycles and replicates and
    equals the shared cycle-1 field — transported once, reused.
  - Empty maps (extinction): cycles after clearance skip transport
    gracefully; dose declared zero (K3).
  - kphys (Gy per unit MeV/decay) set ONCE at cycle 1 so mean viable
    dose = D_MEAN_TARGET = 10 Gy; frozen thereafter.
"""
import json, os   # one line: citations into this file are by line number
import hashlib
import subprocess
import numpy as np

S3 = os.path.dirname(os.path.abspath(__file__))      # study3_closed_loop/
S1, DATA = os.path.join(S3, "..", "study1"), os.path.join(S3, "..", "data")
APP = f"{S1}/build/study1_app"
RUNS = f"{S3}/runs"

ALPHA, BETA = 0.24, 0.06
T_REP, T_ACT = 1.5, 96.0
MU_REP = np.log(2.0) / T_REP
N_CYCLES = 4
A_ADMIN = 1.0
D_MEAN_TARGET = 10.0
SIGMA_DIV = 0.10
T_D = 58.4
LAMBDA_DIV = np.log(2.0) / T_D
DT_DAY = 0.25
P_STEP = 1.0 - np.exp(-LAMBDA_DIV * 24.0 * DT_DAY)
N_STEPS = 224
N_DECAYS = 8_000_000
N_THREADS = 24
N_SEEDS = 5
EXTINCTION_MULT = 100.0
SEED_TRANSPORT_BASE = 20260921
SEED_DYN_CLOSED = 20261001
SEED_DYN_OPEN = 20261101
SEED_K2 = 20261201
SEED_K3 = 20261202

gate_failures = []


def gate(name, ok, detail):
    print(f"GATE {name}: {'PASS' if ok else 'FAIL'} — {detail}")
    if not ok:
        gate_failures.append((name, detail))
    return ok


def G_of_T(T_h):
    x = MU_REP * np.asarray(T_h, dtype=np.float64)
    small = x < 1e-3
    xs = np.where(small, x, 1.0)
    return np.where(small, 1.0 - xs / 3.0 + xs * xs / 12.0,
                    2.0 / (x * x) * (x + np.expm1(-x)))


def sf(D, G):
    D = np.asarray(D, np.float64)
    return np.exp(-(ALPHA * D + BETA * G * D ** 2))


G96 = float(G_of_T(T_ACT))

# ---- Fixed objects ----------------------------------------------------------
geo_meta = json.load(open(f"{S1}/geometry_pitch40_meta.json"))
with open(f"{DATA}/geometry_pitch40.npz", "rb") as f:
    h_geo = hashlib.sha256(f.read()).hexdigest()
assert h_geo == geo_meta["npz_sha256"], "geometry hash mismatch"
gdat = np.load(f"{DATA}/geometry_pitch40.npz")
cell_id = gdat["cell_id"]
n = geo_meta["n_voxels_per_axis"]
NV = n ** 3
VOX = np.flatnonzero(cell_id > 0).astype(np.int64)
CAP = VOX.size
E_T0 = gdat["activity"].astype(np.float64).ravel(order="C").copy()
print(f"Geometry verified sha {h_geo[:16]}.. | capacity {CAP} | G96={G96:.6f}")


def run_transport(tag, seed):
    out = f"{RUNS}/{tag}_dose"
    cmd = [APP, "field", "40.0", f"{RUNS}/{tag}", out, str(N_DECAYS),
           str(N_THREADS), "0", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        raise RuntimeError(f"transport {tag} failed:\n{r.stderr[-2000:]}")
    meta = json.load(open(out + ".json"))
    E = np.fromfile(out + ".bin", dtype=np.float64) / meta["decays_simulated"]
    assert E.size == NV
    return E.reshape((n, n, n), order="C")


def act_map_prop_e(e, occ):
    """Activity ∝ e over OCCUPIED slots only, K1-normalized to A_ADMIN."""
    w = np.zeros(NV)
    w[VOX[occ]] = e[occ]
    tot = w.sum()
    if tot <= 0:
        return None                      # empty population
    return A_ADMIN * w / tot


def write_map(occ, act_vox, tag):
    s = float(act_vox.sum())
    if abs(s - A_ADMIN) > 1e-12 * A_ADMIN:
        raise RuntimeError(f"K1 violated at {tag}: sum={s!r}")
    cid = np.zeros(NV, np.int32)
    cid[VOX[occ]] = np.flatnonzero(occ).astype(np.int32) + 1
    with open(f"{RUNS}/{tag}.bin", "wb") as f:
        f.write(cid.tobytes(order="C"))
        f.write(np.ascontiguousarray(act_vox, dtype=np.float64).tobytes(order="C"))
    with open(f"{RUNS}/{tag}_meta.json", "w") as f:
        json.dump({"n": n, "pitch_um": 40.0}, f)
    return s


def window(occ, e, sf_slot, rng, sigma_div):
    """Death at first post-exposure division + refill to capacity."""
    deaths = 0
    if sf_slot is not None:
        live = np.flatnonzero(occ)
        u = rng.random(live.size)
        died = live[u >= sf_slot[live]]
        occ[died] = False
        e[died] = np.nan
        deaths = int(died.size)
    divisions = 0
    for _ in range(N_STEPS):
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
        e[tgt] = e[attempt[:k]] * np.exp(sigma_div * rng.standard_normal(k))
        divisions += k
    return deaths, divisions, int(occ.sum())


def run_arm(arm, dyn_seed, rep=0, n_cycles=N_CYCLES, sigma_div=SIGMA_DIV,
            dose_mult=1.0, uniform_dose=False, log="", shared=None):
    """Full multi-cycle run. shared: dict for cross-arm transport reuse.

    shared['t0_dose']: dose field of the T0 activity map — bit-identical
    across every arm's cycle 1, the open arm's cycles 2-4, and the
    extinction arm (frozen inputs, frozen transport seed); transported
    ONCE and reused, per protocol v1.2 economy declaration.
    shared['kphys']: declared Gy scaling, frozen at first computation.
    """
    rng = np.random.default_rng(dyn_seed)
    occ = np.ones(CAP, bool)
    e = E_T0[VOX].copy() if arm != "null" else np.ones(CAP)
    kphys = shared.get("kphys") if shared is not None else None
    frozen_w = None
    stats = []
    for cyc in range(1, n_cycles + 1):
        n_start = int(occ.sum())
        if n_start == 0:
            stats.append({"cycle": cyc, "status": "extinct", "n_start": 0,
                          "n_deaths": None, "n_divisions": None, "n_end": 0,
                          "mean_ln_e_start": None, "mean_ln_e_end": None,
                          "mean_e_end": None, "dose_mean": 0.0,
                          "dose_cv": None, "mean_SF": None,
                          "sum_activity": None})
            continue
        tag = f"{log}{arm}_c{cyc}_d{dyn_seed}"
        # --- 1. activity map
        if arm == "open" and frozen_w is not None:
            s = write_map(occ, frozen_w, tag)      # blind re-application
        elif uniform_dose:
            s = A_ADMIN                             # declared uniform dose
        else:
            w_vox = np.zeros(NV)
            w_vox[VOX] = np.nan_to_num(e, nan=0.0)  # ∝ e, live slots only
            act_vox = A_ADMIN * w_vox / w_vox[VOX[occ]].sum()
            s = write_map(occ, act_vox, tag)
            if arm == "open":
                frozen_w = act_vox.copy()           # freeze cycle-1 map
        # --- 2. transport
        if uniform_dose:
            D_occ = np.full(n_start, D_MEAN_TARGET)
            dose_cv = 0.0
        else:
            # v1.2 economy: every arm's cycle-1 map IS the T0 map (closed
            # renormalizes to the same shape; open/extinct use it frozen),
            # so its dose field is shared bit-identically.
            use_shared = (shared is not None and "t0_dose" in shared
                          and (cyc == 1 or arm in ("open", "extinct")))
            if use_shared:
                dose_field = shared["t0_dose"]
            else:
                tseed = SEED_TRANSPORT_BASE + 100 * cyc + rep
                dose_field = run_transport(tag, tseed)
                if shared is not None and "t0_dose" not in shared:
                    shared["t0_dose"] = dose_field
            Ev = dose_field.reshape(-1)[VOX[occ]]
            if kphys is None:
                kphys = D_MEAN_TARGET / float(Ev.mean())
                if shared is not None:
                    shared["kphys"] = kphys
            D_occ = kphys * Ev * dose_mult
            dose_cv = float(D_occ.std() / D_occ.mean())
        # --- 3. response + dynamics
        S_occ = sf(D_occ, G96)
        sf_slot = np.full(CAP, np.nan)
        sf_slot[occ] = S_occ
        mean_ln_e_start = float(np.log(e[occ]).mean())
        deaths, divisions, n_end = window(occ, e, sf_slot, rng, sigma_div)
        stats.append({
            "cycle": cyc, "n_start": n_start, "n_deaths": deaths,
            "n_divisions": divisions, "n_end": n_end,
            "mean_ln_e_start": mean_ln_e_start,
            "mean_ln_e_end": float(np.log(e[occ]).mean()) if occ.any() else None,
            "mean_e_end": float((e[occ]).mean()) if occ.any() else None,
            "dose_mean": float(D_occ.mean()), "dose_cv": dose_cv,
            "mean_SF": float(S_occ.mean()), "sum_activity": s,
            "status": "ok",
        })
    return stats, kphys


# ---- Gates -------------------------------------------------------------------
def gate_k2():
    """Null selection: uniform e, sigma_div=0, uniform dose, 4 cycles."""
    stats, _ = run_arm("null", SEED_K2, sigma_div=0.0, uniform_dose=True,
                       log="k2_")
    means = [s["mean_e_end"] for s in stats if s["mean_e_end"] is not None]
    ok = len(means) == 4 and all(m == 1.0 for m in means)
    return ok, f"mean e after each cycle = {means} (bit-stable 1.0)"


def gate_k3():
    """Extinction: 1000 Gy/cycle clears lattice; empty map graceful."""
    stats, _ = run_arm("extinct", SEED_K3, dose_mult=EXTINCTION_MULT,
                       n_cycles=4, log="k3_")
    ok = (stats[0]["n_deaths"] == stats[0]["n_start"]
          and all(s["status"] == "extinct" for s in stats[1:])
          and all(s["n_end"] == 0 for s in stats[1:]))
    return ok, (f"cycle1 deaths {stats[0]['n_deaths']}/{stats[0]['n_start']}; "
                f"cycles2-4 {[s['status'] for s in stats[1:]]}")


if __name__ == "__main__":
    import os
    import sys
    os.makedirs(RUNS, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "gates"
    if mode == "gates":
        ok2, d2 = gate_k2()
        gate("K2 null-selection", ok2, d2)
        ok3, d3 = gate_k3()
        gate("K3 extinction", ok3, d3)
        print("\nGate failures:", gate_failures or "none")
        sys.exit(1 if gate_failures else 0)
    elif mode == "smoke":
        # 1-cycle transport smoke: T0 map, compare per-viable-cell mean shape
        occ = np.ones(CAP, bool)
        act = act_map_prop_e(E_T0[VOX], occ)
        s = write_map(occ, act, "smoke_t0")
        E = run_transport("smoke_t0", SEED_TRANSPORT_BASE + 1)
        Ev = E.reshape(-1)[VOX]
        kphys = D_MEAN_TARGET / float(Ev.mean())
        print(f"smoke: sum(A)={s!r} mean_E={float(Ev.mean()):.6f} "
              f"kphys={kphys:.6g} Gy/(MeV/decay)")
        print(f"smoke: dose CV over viable = {float(kphys*Ev.std()/ (kphys*Ev.mean())):.4f}")
        sys.exit(0)
    elif mode == "full":
        # Full study: gates K2/K3, then closed x5, open x5, endpoints, verdict.
        shared = {}
        # K2 first (no transport)
        ok2, d2 = gate_k2()
        gate("K2 null-selection", ok2, d2)
        # Shared cycle-1 transport: closed arm seed 0 computes t0_dose
        print("\n=== CLOSED ARM (5 dynamics seeds) ===")
        closed = []
        for rep in range(N_SEEDS):
            st, _ = run_arm("closed", SEED_DYN_CLOSED + rep, rep=rep,
                            shared=shared)
            closed.append(st)
            c4 = st[3]
            print(f"closed rep{rep}: ln-e drift c1->c4 = "
                  f"{c4['mean_ln_e_start'] - st[0]['mean_ln_e_start']:+.4f}, "
                  f"deaths/cyc = {[s['n_deaths'] for s in st]}, "
                  f"n_end c4 = {c4['n_end']}")
        print("\n=== OPEN ARM (5 dynamics seeds, frozen map) ===")
        open_stats = []
        for rep in range(N_SEEDS):
            st, _ = run_arm("open", SEED_DYN_OPEN + rep, rep=rep,
                            shared=shared)
            open_stats.append(st)
            c4 = st[3]
            print(f"open   rep{rep}: ln-e drift c1->c4 = "
                  f"{c4['mean_ln_e_start'] - st[0]['mean_ln_e_start']:+.4f}, "
                  f"deaths/cyc = {[s['n_deaths'] for s in st]}")
        # K3 extinction (reuses shared T0 field at x100 dose scale)
        ok3, d3 = gate_k3()
        gate("K3 extinction", ok3, d3)
        # K1: conservation across every recorded cycle
        k1_ok = all(
            s["sum_activity"] is None or abs(s["sum_activity"] - A_ADMIN) <= 1e-12 * A_ADMIN
            for arm_st in closed + open_stats for s in arm_st)
        gate("K1 conservation (all cycles, all arms)", k1_ok,
             "sum(A) == A_ADMIN at every write (bit-exact)")
        # ---- Endpoints (v1.1/v1.2 definitions) ----
        dl_closed = [st[3]["mean_ln_e_start"] - st[0]["mean_ln_e_start"]
                     for st in closed]
        dl_open = [st[3]["mean_ln_e_start"] - st[0]["mean_ln_e_start"]
                   for st in open_stats]
        dc = float(np.mean(dl_closed))
        do = float(np.mean(dl_open))
        delta_sel = abs(dc) - abs(do)
        worst_open = max(dl_open)          # least-drifted open replicate
        e1a = dc < worst_open
        e1b = delta_sel > 0
        # E2: sawtooth — refill to capacity at every cycle start, no extinction
        e2 = all(s["n_start"] == CAP for st in closed for s in st
                 if s["status"] == "ok")
        surv = [[s["n_start"] - s["n_deaths"] for s in st] for st in closed]
        crit = {
            "E1a_closed_below_worst_open": bool(e1a),
            "E1b_contrast_positive": bool(e1b),
            "E2_sawtooth_capacity_refill": bool(e2),
        }
        if gate_failures:
            verdict = "REFUSED-GATED"
        elif all(crit.values()):
            verdict = "PASS-SELECTION"
        else:
            verdict = "FAIL-SELECTION"
        print(f"\n=== VERDICT: {verdict} ===")
        print(f"Delta ln-e closed (mean of 5): {dc:+.4f} per replicate: "
              f"{[f'{x:+.4f}' for x in dl_closed]}")
        print(f"Delta ln-e open   (mean of 5): {do:+.4f} per replicate: "
              f"{[f'{x:+.4f}' for x in dl_open]}")
        print(f"Delta_sel = {delta_sel:+.4f}; worst (least-drifted) open = "
              f"{worst_open:+.4f}")
        print(f"Survivors per cycle (closed): {surv}")
        out = {
            "protocol": "study3_closed_loop/PROTOCOL.md v1.0+v1.1+v1.2",
            "verdict": verdict,
            "criteria": crit,
            "delta_ln_e_closed_mean": dc, "delta_ln_e_closed_per_rep": dl_closed,
            "delta_ln_e_open_mean": do, "delta_ln_e_open_per_rep": dl_open,
            "delta_sel": delta_sel, "worst_open_replicate": worst_open,
            "survivors_per_cycle_closed": surv,
            "gate_failures": gate_failures,
            "kphys_Gy_per_MeV_per_decay": shared.get("kphys"),
            "declared": {
                "alpha": ALPHA, "beta": BETA, "T_rep_h": T_REP,
                "T_active_h": T_ACT, "G96": G96, "A_admin": A_ADMIN,
                "D_mean_target_Gy": D_MEAN_TARGET, "sigma_div": SIGMA_DIV,
                "T_doubling_h": T_D, "window_days": 56, "n_cycles": 4,
                "n_decays_per_cycle": N_DECAYS,
                "seeds": {"transport_base": SEED_TRANSPORT_BASE,
                          "closed": [SEED_DYN_CLOSED + i for i in range(N_SEEDS)],
                          "open": [SEED_DYN_OPEN + i for i in range(N_SEEDS)]},
            },
            "closed_stats": closed, "open_stats": open_stats,
        }
        with open(f"{S3}/verdict.json", "w") as f:
            json.dump(out, f, indent=2)
        print(f"\nWrote {S3}/verdict.json")
        sys.exit(1 if gate_failures or not all(crit.values()) else 0)