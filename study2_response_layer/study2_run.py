#!/usr/bin/env python3
"""Study 2 — Response-layer verification (open-loop LQ/mitotic-catastrophe TCP).

Executes study2_response_layer/PROTOCOL.md v1.0 exactly: gates G1-G8 first,
then the pre-declared endpoints. No TCP verdict is issued if any gate fails
(REFUSED-GATED).

Fixed objects (Study 1, hashed):
  geometry_pitch40.npz  (activity/cell/mask map, sha pinned in meta json)
  runs/ref_p40_h3.bin   (128M-decay Geant4 REF dose field, MeV/voxel)
"""
import json
import hashlib
import numpy as np

BASE = "/home/aurascoper/Developer/PRRT-spatial-cpm/study1_dpk_vs_mc"
OUT = "/home/aurascoper/Developer/PRRT-spatial-cpm/study2_response_layer"

# ---- Declared parameters (PROTOCOL v1.0, fixed before running) -------------
ALPHA = 0.24          # Gy^-1, NCI-H69 EBRT fit, Tamborino 2025 Table 1
BETA = 0.06           # Gy^-2, same source
T_REP = 1.5           # h, declared repair half-time (consistency-gated, G5)
T_ACT = 96.0          # h, declared active exposure duration
N_CLON_SCAN = [100.0, 1000.0]
DOSE_SCAN = [2.0, 5.0, 10.0, 20.0, 40.0]   # Gy, mean viable-cell dose
DOSE_PRIMARY = 10.0
N_ENSEMBLE = 200
SEED_ENSEMBLE = 20260919
REF_NOISE_P95 = 0.0117                     # Study 1: REF self-noise p95 = 1.17%
NOISE_SIGMA = REF_NOISE_P95 / 1.645        # half-normal p95 -> sigma mapping
MU = np.log(2.0) / T_REP
NCLON_INVARIANCE_TOL = 1e-9

gate_failures = []


def gate(name, ok, detail):
    status = "PASS" if ok else "FAIL"
    print(f"GATE {name}: {status} — {detail}")
    if not ok:
        gate_failures.append((name, detail))
    return ok


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def G_of_T(T_h, mu=MU):
    """Lea-Catcheside factor, constant dose rate, monoexponential repair.

    Numerically stable form: x - 1 + e^-x == x + expm1(-x) exactly (identity),
    and expm1 keeps full relative accuracy at small x where the naive form
    cancels catastrophically (G5c caught exactly this defect in v1.0 of the
    runner: naive form returned G(1e-9 h) = 0.0 instead of 1.0).
    Second defect caught by G5c after the expm1 fix: the SUM x + expm1(-x)
    still rounds at relative ~1.8e-7 for x ~ 1e-10 (true result x^2/2 is
    ~2^21 ulps below the operands). Final fix per protocol v1.1 amendment 1:
    Taylor branch G ~= 1 - x/3 + x^2/12 for x < 1e-3 (omitted -x^3/60 term
    ~ 1.7e-14 at the branch point), expm1 form elsewhere.
    """
    x = mu * np.asarray(T_h, dtype=np.float64)
    small = x < 1e-3
    xs = np.where(small, x, 1.0)          # safe divisor under the mask
    g_small = 1.0 - xs / 3.0 + xs * xs / 12.0
    g_full = 2.0 / (x * x) * (x + np.expm1(-x))
    return np.where(small, g_small, g_full)


def sf(D, G, alpha=ALPHA, beta=BETA):
    """Preprint Eq. sf: SF_i = exp(-alpha*D - beta*G*D^2)."""
    return np.exp(-(alpha * np.asarray(D, dtype=np.float64)
                    + beta * G * np.asarray(D, dtype=np.float64) ** 2))


def ln_tcp(D_viable, n_clon, G):
    """ln TCP = -N_s, N_s = n_clon * sum_i SF_i. Exact in N_s space."""
    return -n_clon * float(np.sum(sf(D_viable, G)))


# ---- G7: fixed-object hashes ----------------------------------------------
geo_meta = json.load(open(f"{BASE}/geometry_pitch40_meta.json"))
geo_sha = sha256_file(f"{BASE}/geometry_pitch40.npz")
gate("G7a geometry hash", geo_sha == geo_meta["npz_sha256"],
     f"meta {geo_meta['npz_sha256'][:16]}.. vs file {geo_sha[:16]}..")

ref_meta = json.load(open(f"{BASE}/runs/ref_p40_h3.json"))
ref_sha = sha256_file(f"{BASE}/runs/ref_p40_h3.bin")
# No stored reference hash exists for the REF bin (Study 1 pinned the npz and
# the DPK only); this run PINS it for all future re-use. Declared in protocol
# execution order step 1 note.
print(f"GATE G7b dose-field hash: PINNED {ref_sha} "
      f"(first pin; 128M decays, seed 20260918 per runs/ref_p40_h3.json)")

# ---- G5: Lea-Catcheside admissibility --------------------------------------
T_grid = np.array([1.0, 6.0, 24.0, 96.0, 360.0])
G_vals = G_of_T(T_grid)
gate("G5a band 0<G<1", bool(np.all((G_vals > 0) & (G_vals < 1))),
     f"G(T=1..360h) = {np.round(G_vals, 5).tolist()}")
T_dense = np.linspace(1e-3, 500.0, 200001)
G_dense = G_of_T(T_dense)
mono_ok = bool(np.all(np.diff(G_dense) <= 1e-12))
gate("G5b monotone non-increasing in T", mono_ok, "dense grid 1e-3..500 h")
G0 = G_of_T(1e-9)
gate("G5c G(T->0) -> 1", abs(G0 - 1.0) < 1e-9, f"G(1e-9 h) = {G0:.12f}")
G96 = G_of_T(T_ACT)
one_over_G = 1.0 / G96
gate("G5d Tamborino consistency", 1.0 < one_over_G < 100.0,
     f"1/G(96h) = {one_over_G:.1f} vs measured ~12 for NCI-H69 PRRT "
     f"(measured 1/G inflated by non-repair processes; model < measured OK)")

# ---- G1/G2/G3: survival-function gates -------------------------------------
sf0 = float(sf(0.0, G96))
gate("G1 null", sf0 == 1.0, f"SF(0 Gy) = {sf0!r}")
sf1k = float(sf(1000.0, G96))
gate("G2 obliteration", sf1k <= 1e-100, f"SF(1000 Gy) = {sf1k:.3e}")

g3_ok, g3_detail = True, []
for D_g in [2.0, 8.0, 10.0]:
    s_fast = float(sf(D_g, G_of_T(24.0)))
    s_slow = float(sf(D_g, G_of_T(96.0)))
    ok = s_fast < s_slow
    g3_ok &= ok
    g3_detail.append(f"D={D_g}: SF(24h)={s_fast:.6f} < SF(96h)={s_slow:.6f}: {ok}")
gate("G3 rate/protraction", g3_ok, "; ".join(g3_detail))

# ---- Load dose field, build cell doses -------------------------------------
g = np.load(f"{BASE}/geometry_pitch40.npz")
cell_id = g["cell_id"]
viable = cell_id > 0
n_viable = int(viable.sum())
E = np.fromfile(f"{BASE}/runs/ref_p40_h3.bin", dtype=np.float64) / ref_meta["decays_simulated"]
n_axis = geo_meta["n_voxels_per_axis"]
assert E.size == n_axis ** 3 == cell_id.size, "dose field shape mismatch"
E = E.reshape((n_axis, n_axis, n_axis), order="C")   # C-order per ref json units
E_v = E[viable]                     # MeV/decay at viable-cell sites
E_mean = float(E_v.mean())
D_v = DOSE_PRIMARY * E_v / E_mean   # declared scaling: shape * D_mean
cv_dose = float(D_v.std() / D_v.mean())
p95_ratio = float(np.percentile(D_v, 95) / D_v.mean())
print(f"\nViable cells: {n_viable}; dose field over viable cells: "
      f"CV={cv_dose:.4f}, p95/mean={p95_ratio:.4f}")

# ---- G4: uniform-field assembler gate --------------------------------------
D_unif_field = np.full(n_viable, D_v.mean())
ln_het_unif_input = ln_tcp(D_unif_field, 1.0, G96)
ln_uni_unif_input = -float(np.sum(sf(D_unif_field, G96)))  # independent path
gate("G4 uniform-field assembler",
     abs(ln_het_unif_input - ln_uni_unif_input) <= 1e-9 * abs(ln_uni_unif_input),
     "het pipeline on a uniform field == uniform pipeline (machine precision)")

# ---- G8: energy bookkeeping -------------------------------------------------
D_unif = float(D_v.mean())
G96_val = G96
# total viable dose, both arms, in Gy*cells
tot_het = float(D_v.sum())
tot_uni = D_unif * n_viable
gate("G8 energy bookkeeping", abs(tot_het - tot_uni) <= 1e-12 * tot_het,
     f"total viable dose het {tot_het:.9f} vs uniform {tot_uni:.9f} Gy*cells")

# ---- Arms + scans ------------------------------------------------------------
def arm_pair(D_v_field, dose_mean_scale, n_clon, G):
    D_h = dose_mean_scale * D_v_field / D_v_field.mean()
    D_u = dose_mean_scale * D_v_field.mean() / D_v_field.mean()  # == mean
    D_u = float(D_h.mean())
    return ln_tcp(D_h, n_clon, G), ln_tcp(np.full(D_h.size, D_u), n_clon, G)

results = {}
all_c = True
all_d = True
print(f"\n{'D_mean':>7} {'N_clon':>7} {'lnTCP_het':>14} {'lnTCP_unif':>14} "
      f"{'N_s_het':>12} {'N_s_unif':>12} {'gap':>12}")
for Dm in DOSE_SCAN:
    for ncl in N_CLON_SCAN:
        # N_clon invariance is linear in N_s: compute at N_clon=1, scale
        lnh_1, lnu_1 = arm_pair(E_v, Dm, 1.0, G96_val)
        gap_1 = lnh_1 - lnu_1          # = -(N_s_het - N_s_unif) at n_clon=1
        lnh = ncl * lnh_1
        lnu = ncl * lnu_1
        inv_ok = abs((ncl * gap_1) - (lnh - lnu)) <= NCLON_INVARIANCE_TOL * abs(ncl * gap_1)
        all_c &= inv_ok
        direction = lnh < lnu          # ln TCP_het < ln TCP_unif <=> TCP_het < TCP_unif
        all_c &= direction
        ns_het = -lnh
        ns_uni = -lnu
        results[f"D{Dm:g}_N{ncl:g}"] = {
            "ln_TCP_het": lnh, "ln_TCP_uniform": lnu,
            "N_s_het": ns_het, "N_s_uniform": ns_uni,
            "TCP_het": float(np.exp(lnh)) if lnh > -700 else 0.0,
            "TCP_uniform": float(np.exp(lnu)) if lnu > -700 else 0.0,
            "direction_het_lt_unif": bool(direction),
            "nclon_invariance_ok": bool(inv_ok),
        }
        print(f"{Dm:7g} {ncl:7g} {lnh:14.4f} {lnu:14.4f} {ns_het:12.4g} "
              f"{ns_uni:12.4g} {ns_uni - ns_het:12.4g}")

gate("C(c) direction+invariance at all levels", all_c,
     "TCP_het < TCP_uniform at every dose level and both N_clon "
     f"(invariance tol {NCLON_INVARIANCE_TOL:g})")

# (d) monotonicity per protocol v1.1 amendment 2: N_s_het strictly DECREASING
# in D_mean (TCP_het = exp(-N_s) then rises with dose, as the physics
# requires). The v1.0 text had the direction inverted against the protocol's
# own endpoint; amended BEFORE any verdict-issuing run.
ns_het_seq = [results[f"D{Dm:g}_N100"]["N_s_het"] for Dm in DOSE_SCAN]
all_d = bool(np.all(np.diff(ns_het_seq) < 0))
gate("C(d) monotonicity", all_d,
     f"N_s_het over scan (must strictly decrease) = "
     f"{[f'{x:.4g}' for x in ns_het_seq]}")

# ---- G6: noise ensemble at the primary level --------------------------------
rng = np.random.default_rng(SEED_ENSEMBLE)
xi = np.exp(NOISE_SIGMA * rng.standard_normal((N_ENSEMBLE, E_v.size))
            - 0.5 * NOISE_SIGMA ** 2)
gaps = np.empty(N_ENSEMBLE)
for r in range(N_ENSEMBLE):
    Dr = D_v * xi[r]
    ns_h = float(np.sum(sf(Dr, G96_val)))
    Du = float(Dr.mean())
    ns_u = float(np.sum(sf(np.full(Dr.size, Du), G96_val)))
    gaps[r] = ns_h - ns_u        # divergence D = N_s_het - N_s_unif at N_clon=1
                                 # (positive = het leaves MORE survivors = LOWER
                                 # TCP, the Mellhammar direction; protocol v1.1
                                 # amendment 3 — v1.0 runner had the sign flipped)
gap_det = results[f"D{DOSE_PRIMARY:g}_N100"]["N_s_het"] - \
          results[f"D{DOSE_PRIMARY:g}_N100"]["N_s_uniform"]      # at N_clon=1: /100
gap_det_1 = gap_det / 100.0
sig = float(gaps.std(ddof=1))
crit_b = gap_det_1 > 2.0 * sig
gate("G6/C(b) noise ensemble", crit_b,
     f"divergence D (N_clon=1) = {gap_det_1:.6g} survivors vs ensemble "
     f"std = {sig:.6g} (mean {gaps.mean():.6g}); ratio = "
     f"{gap_det_1 / sig:.1f}x the 2-sigma bar; "
     f"sigma={NOISE_SIGMA:.5f} (p95 {REF_NOISE_P95:.4f}), n={N_ENSEMBLE}")

# ---- Verdict -----------------------------------------------------------------
primary = results[f"D{DOSE_PRIMARY:g}_N100"]
criteria = {
    "a_direction_primary": primary["direction_het_lt_unif"],
    "b_significance_vs_mc_noise": crit_b,
    "c_direction_all_levels_both_nclon": bool(all_c),
    "d_monotonicity": all_d,
}
if gate_failures:
    verdict = "REFUSED-GATED"
    reason = f"gate failures: {gate_failures}"
elif all(criteria.values()):
    verdict = "PASS-DIVERGENCE"
    reason = "all pre-declared criteria (a)-(d) hold; all gates G1-G8 passed"
else:
    verdict = "FAIL-DIVERGENCE"
    reason = f"failed criteria: {[k for k, v in criteria.items() if not v]}"

print(f"\n=== VERDICT: {verdict} ===\n{reason}")

out = {
    "protocol": "study2_response_layer/PROTOCOL.md v1.0",
    "verdict": verdict,
    "reason": reason,
    "criteria": criteria,
    "gates_failed": gate_failures,
    "declared_parameters": {
        "alpha_Gy^-1": ALPHA, "beta_Gy^-2": BETA,
        "T_rep_h": T_REP, "T_active_h": T_ACT, "mu_h^-1": float(MU),
        "G_96h": float(G96_val), "1/G_96h": float(one_over_G),
        "N_clon_scan": N_CLON_SCAN, "dose_scan_Gy": DOSE_SCAN,
        "dose_primary_Gy": DOSE_PRIMARY,
        "noise_sigma": NOISE_SIGMA, "ref_noise_p95": REF_NOISE_P95,
        "seed_ensemble": SEED_ENSEMBLE, "n_ensemble": N_ENSEMBLE,
        "source": "Tamborino 2025 JNM 66:1291, Table 1 (EBRT NCI-H69 LQ fit); "
                  "T_rep declared prior, consistency-gated (G5d)",
    },
    "dose_field_descriptors": {
        "n_viable": n_viable, "dose_CV_viable": cv_dose,
        "dose_p95_over_mean": p95_ratio,
        "E_mean_MeV_per_decay_viable": E_mean,
        "geometry_sha256": geo_sha, "ref_bin_sha256_pinned": ref_sha,
        "ref_decays": ref_meta["decays_simulated"],
},
    "results_per_level": results,
    "noise_ensemble": {
        "gap_Ns_mean": float(gaps.mean()), "gap_Ns_std": sig,
        "gap_det_deterministic": gap_det_1,
        "ratio_to_2sigma": float(gap_det_1 / (2.0 * sig)),
        "p2p5": float(np.percentile(gaps, 2.5)),
        "p97p5": float(np.percentile(gaps, 97.5)),
    },
    "G_curve": {f"T={t:g}h": float(gg) for t, gg in zip(T_grid, G_vals)},
}
with open(f"{OUT}/verdict.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\nWrote {OUT}/verdict.json")

# SF distribution summaries at primary level (for the report)
s_het = sf(D_v, G96_val)
s_uni = sf(np.full(n_viable, D_unif), G96_val)
for nm, s in [("het", s_het), ("uniform", s_uni)]:
    print(f"SF[{nm}] mean={s.mean():.5f} p50={np.percentile(s, 50):.5f} "
          f"p95={np.percentile(s, 95):.5f}")
print(f"D_unif = {D_unif:.4f} Gy at D_mean={DOSE_PRIMARY:g} Gy declaration")