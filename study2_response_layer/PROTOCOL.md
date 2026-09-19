# Response-Layer Verification Study — Protocol

Project: PRRT-spatial-cpm, Step 2 of the dependency order.
Status: PROTOCOL (pre-declared, before any TCP number is computed).
Written: 2026-09-17. This document precedes all results; every parameter,
gate, and verdict criterion below was fixed before the first survival
evaluation was run. Build on Study 1 (FAIL-DPK verdict, protocol v1.5,
commits 2319b7a/d715000) as a trusted transport foundation.

## Purpose

The closed-loop PRRT model (preprint v1.1, Sec. "Response: protracted-exposure
survival kinetics") must reproduce, in open loop, the tumor-control
divergence that Mellhammar et al. established (JNM 2023;64:1632): a static
heterogeneous intratumoral activity distribution yields a different TCP than
an assumed uniform dose at equal total energy. This study verifies that the
declared response layer — protracted Linear-Quadratic survival (Eqs. sf/lea)
over the validated cell-resolved dose field — mathematically reproduces that
divergence, with pre-declared gates that make a broken implementation fail
loudly.

What this study is NOT: it is not a TCP validation against clinical data.
Both arms are model constructs on a synthetic geometry. The deliverable is
the verified existence, direction, and noise-robustness of the divergence,
on the exact transport field Study 1 validated.

## Fixed objects (hashed; consumed, never regenerated)

1. Geometry/activity: study1_dpk_vs_mc/geometry_pitch40.npz,
   sha256 = a6883b84a577dd5db888c9d4292cf68897a68d840b64a427061bc70a2d60a170
   (from geometry_pitch40_meta.json). 25^3 voxels at 40 um pitch, 13,978
   viable cells (cell_id > 0), tissue masks, lognormal activity map
   (sigma_e = 0.6, activity CV = 0.658). Same hashed object as Study 1.
2. Dose field: study1_dpk_vs_mc/runs/ref_p40_h3.bin (128M decays, seed
   20260918, Geant4 11.4.2, ENSDF RDM 6.1.2), units MeV/voxel, C-order,
   per-decay normalization by 128e6. REF self-noise p95 = 1.17% at the cell
   endpoint (Study 1, protocol v1.5).
3. The DPK convolution arm is EXCLUDED from this study by the Study-1
   FAIL-DPK verdict. The response layer runs only on REF re-transport, as
   mandated by protocol v1.5 consequence 1.

A hash mismatch on either object aborts the run.

## Declared radiobiological parameters (fixed before running)

Source: Tamborino et al., JNM 2025;66:1291-1298
(doi:10.2967/jnumed.125.269470), Table 1 (best-fitting dose-response
parameters) and Discussion (1/G analysis).

1. Intrinsic LQ coefficients: the ACUTE (EBRT) fit for the proliferative
   NET line NCI-H69 — alpha = 0.24 Gy^-1 (0.24 +- 0.05), beta = 0.06 Gy^-2
   (0.06 +- 0.02), alpha/beta = 4 Gy.
   [Rationale, declared: the PRRT-fitted coefficients of the same paper
   (alpha = 0.03 +- 0.10, beta = 0.09 +- 0.07) ALREADY embed the protraction
   effect through their fitted dose rate. Applying our separate Lea-
   Catcheside factor on top of them would double-count repair. The EBRT-fit
   parameters are the intrinsic radiosensitivity to which G(T) is applied.
   This is the standard decomposition and is consistent with Tamborino's
   own treatment, which quantifies protraction via 1/G.]
2. Exposure protraction: constant dose rate over active duration T = 96 h
   (4 days; Tamborino reports most of the absorbed dose delivered within
   3-4 days, effective half-life 22.6 +- 0.15 h for NCI-H69; our synthetic
   exposure uses declared constant-rate idealization, and the deviation is
   noted as a declared modelling choice, class D).
3. Repair half-time: T_rep = 1.5 h (mu = ln 2 / T_rep).
   [Declared literature prior for sublethal-damage repair in the LQ
   framework; Tamborino 2025 publishes no explicit T_rep. The protocol
   therefore declares a CONSISTENCY GATE rather than silent adoption: with
   T = 96 h the model's 1/G must be within one order of magnitude of the
   ~12 measured for NCI-H69 PRRT (their 1/G is inflated above pure repair
   by redistribution/reoxygenation/proliferation, so 1/G_model < 12 is
   acceptable).]
4. Lea-Catcheside factor (preprint Eq. lea), constant dose rate:
   G(T) = 2/(mu T)^2 * (mu T - 1 + e^(-mu T)).
   Declared admissibility band: 0 < G < 1, G -> 1 as T -> 0, G monotone
   non-increasing in T (verified numerically before any TCP is computed).
5. Survival (preprint Eq. sf): SF_i = exp(-alpha*D_i - beta*G*D_i^2),
   D_i in Gy. Death is realized at mitosis entry (mitotic catastrophe);
   in this open-loop study the population-level consequence is the
   Poisson TCP below — the stochastic CPM realization is Step 3.
6. Clonogen burden: N_clon = 100 clonogens per viable cell (declared
   choice, class D). The protocol pre-declares the check that the TCP
   divergence DIRECTION is invariant to N_clon over the scan {1e2, 1e3};
   N_clon changes TCP steepness, not the sign of the comparison.
7. Dose scale: the dose field is a shape (MeV/voxel at declared mean
   activity 1 Bq/voxel). One declared, non-fitted conversion constant sets
   the mean viable-cell dose D_mean; the study is run at TWO pre-declared
   levels: D_mean = 10 Gy (primary) and a scan {2, 5, 10, 20, 40} Gy.
   Divergence direction and significance are claimed at every level,
   not a single one. Conversion: D_i (Gy) = E_i (MeV/decay) * A_eff (Bq) *
   T_eff (s) * 1.602176634e-13 (J/MeV) / m_voxel (kg), m_voxel = (40 um)^3
   * 1 g/cm^3 = 6.4e-11 kg, with the time-integration constant folded into
   the declared D_mean normalization (the shape is what carries the
   physics; the absolute scaling is declared, not inferred).
   Cell dose convention retained from Study 1: cell = voxel site at pitch
   40 um; the manuscript states cell-dose = voxel-dose here.

## Arms

- TCP_het: SF evaluated per cell from the heterogeneous REF dose field.
  N_s = N_clon * sum_i SF_i over viable cells; TCP = exp(-N_s).
- TCP_uniform: SF evaluated from a FLAT dose field assigned to every
  viable cell, calibrated to deliver the EXACT same total absorbed energy
  to the viable population as the heterogeneous field:
  D_unif = sum_i D_i (viable) / n_viable. Both arms thus integrate the
  same total dose over the same cell population — the comparison isolates
  the effect of heterogeneity alone.

## PASS/FAIL CRITERIA (pre-declared)

Primary endpoint: TCP_het vs TCP_uniform at D_mean = 10 Gy.

PASS-DIVERGENCE: ALL of
  (a) TCP_het < TCP_uniform at D_mean = 10 Gy (direction fixed by
      pre-declaration: hot cells waste dose on already-doomed clonogens
      through the quadratic term while cold cells survive; Mellhammar's
      TCP dropped under heterogeneous activity distribution);
  (b) the gap is significant against the transport noise floor: with the
      dose-field Monte Carlo noise ensemble defined below, the gap must
      exceed 2 ensemble standard deviations of the gap;
  (c) direction holds at EVERY level of the dose scan {2, 5, 10, 20, 40} Gy
      and for BOTH N_clon values;
  (d) TCP_het is monotone non-increasing in D_mean over the scan.

FAIL-DIVERGENCE: any of (a)-(d) fails; the failure is reported as-is with
the failing arm and level, and the response layer is not declared verified.

REFUSE: any gate in the next section fails; no TCP verdict is issued from a
gated run (refuse, do not pass — same rule as Study 1 C1).

## Negative controls and gates (every check can fail)

G1 — Null gate: a cell receiving D = 0.0 Gy must evaluate to
     SF = exp(0) = 1.0 EXACTLY (bit-exact 1.0 after the exponential).
G2 — Obliteration gate: a cell receiving D = 1000 Gy must evaluate to
     SF <= 1e-100 (computationally indistinguishable from 0; guards
     against exp() underflow-to-one sign/logic errors).
G3 — Rate/protraction gate (Lea-Catcheside): for a FIXED total dose D,
     reducing T (higher dose rate) MUST yield lower SF than extending T.
     Verified numerically at (D, T) grid {(2, 24), (2, 96), (8, 24),
     (8, 96), (10, 96)} Gy-h: SF(D, T=24h) < SF(D, T=96h) must hold at
     every grid point (G(24h) > G(96h) with T_rep = 1.5 h). A broken
     G(T) (sign error, G > 1, wrong limit) fails here.
G4 — Uniform-field gate: when the heterogeneous field is replaced by a
     uniform field at the same mean (the arms collapse to the same input),
     TCP_het and TCP_uniform must agree to machine precision. Guards the
     TCP assembler against per-cell summation defects.
G5 — G(T) admissibility band: 0 < G < 1 over T in {1, 6, 24, 96, 360} h;
     G(T->0) -> 1 within 1e-9; monotone non-increasing over a dense T
     grid. Also the Tamborino consistency check: 1/G(96 h) must lie in
     (1, 100) given the measured 1/G ~ 12 with non-repair contributors
     inflating it (declared order-of-magnitude gate).
G6 — Noise ensemble (MC-noise negative control): resample the cell-dose
     field 200 times with multiplicative lognormal noise whose p95 equals
     the Study-1 REF self-noise (p95 = 1.17% -> sigma = p95/1.645 under a
     half-normal mapping); recompute the gap per resample. The verdict's
     significance test (criterion b) REQUIRES this ensemble; if the gap
     under noise overlaps zero at the 2-sigma level, the study reports
     FAIL-DIVERGENCE (criterion b failed), not a silent pass.
G7 — Hash equality on both fixed objects (abort on mismatch).
G8 — Energy bookkeeping: the uniform arm's total viable-population dose
     must equal the heterogeneous arm's to <= 1e-12 relative (it is
     constructed, not transported; any discrepancy is an implementation
     bug).

## Outputs

- verdict (PASS-DIVERGENCE / FAIL-DIVERGENCE / REFUSED-GATED) with the
  failing gate named;
- TCP_het, TCP_uniform, gap, and the 2-sigma noise band, per dose level
  and per N_clon;
- SF distribution summary (mean, p50, p95) per arm;
- G(T) curve values at the declared T grid;
- dose-shape descriptors of the heterogeneous field over viable cells
  (mean, CV, p95/mean) to tie the divergence magnitude to the measured
  heterogeneity;
- all gate results with measured values.

## Boundaries

- No CPM dynamics, no proliferation, no reoxygenation, no phase-dependent
  radiosensitivity: single-cycle open-loop response only (Step 3 adds the
  closed loop).
- No oxygenation modulation of alpha (no oxygen map in the synthetic
  geometry; preprint defers this to measured-input lesions).
- Phase-dependence of alpha_phi/beta_phi omitted: the declared prior is
  the population-averaged EBRT fit. Declared as a limitation, not a
  silent simplification.
- Constant-dose-rate idealization of a decaying source: T = 96 h active
  duration declared; using T_eff-based effective exposure rather than a
  full biexponential integration is a declared class-D choice.
- TCP = exp(-N_s) is the Poisson clonogen endpoint; it is the endpoint
  Mellhammar-style comparisons use, and the manuscript states that the
  stochastic mitotic-catastrophe realization lives in the CPM (Step 3).

## Execution order (all criteria fixed above, before results)

1. Verify hashes of both fixed objects (G7).
2. Verify G(T) band + monotonicity + Tamborino consistency (G5).
3. Verify functional gates G1-G3 on the survival function alone.
4. Load dose field, build cell doses, both arms; verify G4, G8.
5. Run noise ensemble (G6) and the dose/N_clon scans; evaluate criteria
   (a)-(d).
6. Write verdict.json + this protocol co-committed with the code.

## v1.1 amendments (2026-09-17, issued after the first gated run REFUSED,
## before any verdict-issuing run; no TCP verdict was interpreted from the
## refused run)

Three defects were caught by the protocol's own gates — the discipline
worked; none of them is the physics. Recorded in the Study-1 style:

1. G5c numerical defect, fixed: the Lea-Catcheside implementation
   `2/(x^2) * (x - 1 + e^-x)` cancels catastrophically at small x
   (x = mu*T = 4.6e-10 at the declared T = 1e-9 h probe): even with
   expm1, the SUM x + expm1(-x) rounds at relative ~1.8e-7 because the
   true result x^2/2 is ~2^21 ulps below the operands. First run returned
   G(1e-9 h) = 1.000000176665 and the gate REFUSED the run. Fix: Taylor
   branch G(x) ~= 1 - x/3 + x^2/12 for x < 1e-3 (next omitted term
   -x^3/60 ~ 1.7e-14 at the branch point, far inside the 1e-9 gate),
   expm1 form elsewhere. Numerical-implementation fix, not a threshold
   change; gate and probe unchanged.
2. C(d) direction defect in THIS protocol's own text: v1.0 declared
   "TCP_het monotone non-increasing in D_mean", which contradicts the
   protocol's own endpoint TCP = exp(-N_s): more dose lowers SF_i, hence
   N_s falls and TCP rises. The first run measured N_s_het strictly
   DECREASING over the scan (8.6e5 -> 356) and the gate as written
   refused a physically correct field. Amendment: the pre-declared
   monotonicity check is N_s_het strictly decreasing in D_mean
   (equivalently TCP_het non-decreasing). Endpoint-correcting, not
   threshold-tuning; the direction constraint (a) is untouched.
3. G6 sign-convention defect in the runner: the ensemble gap was computed
   as N_s_uniform - N_s_het, but per criterion (a) the divergence is
   N_s_het > N_s_uniform (heterogeneous field leaves MORE surviving
   clonogens, hence lower TCP). Criterion (b) is evaluated on the
   divergence magnitude D = N_s_het - N_s_uniform > 2 ensemble sigmas.
   The refused run already showed D = 249.36 (N_clon=1 units) against
   ensemble std 0.107 — the sign bug masked a 2330x margin. Sign
   correction only; criterion (b) unchanged.
4. NOTE recorded for the manuscript (from the refused run's numbers, no
   verdict implied): at the declared clonogen burdens, TCP itself
   underflows to 0 for both arms at D_mean <= 20 Gy (N_s ~ 1e4-1e6); the
   divergence is resolved in log-TCP (N_s) space, and in direct TCP space
   only at the high-dose end of the scan (40 Gy: TCP_het ~ e^-356 vs
   TCP_unif ~ 0.25). The Jensen mechanism is visible in the SF summaries:
   E[SF] over the heterogeneous field (0.0875) exceeds SF at the same
   mean dose (0.0696) at 10 Gy because SF(D) is convex there
   ((alpha + 2*beta*G*D)^2 > 2*beta*G), so heterogeneity RAISES mean
   clonogen survival and LOWERS TCP relative to the uniform assumption —
   the Mellhammar direction, now with its mechanism named.
## v1.2 amendment (2026-09-18, after the verdict; issues #11 and #13)

What the triage of 2026-09-18 found. Criteria (a), (c) and (d) follow from the
LQ form for any dose field that is not exactly constant: SF is convex for all
D >= 0 because alpha^2 = 0.0576 exceeds 2*beta*G = 0.00529, so Jensen's
inequality fixes their direction. Criterion (b) compares the gap with its own
spread under re-noising, which measures the gap's precision. Substituting the
REF field with a flat field carrying 0.71% Monte Carlo noise gave
PASS-DIVERGENCE at 48.7x criterion (b); random garbage gave 7503x. Gate G7b
printed the hash it found and compared it with nothing. Gate G8 compared a sum
with mean times count, an identity.

What changes. Three gates, no criterion:

- G6b, noise-only null: the same 200-draw ensemble applied to a flat field at
  the same mean dose. The real gap must exceed twice the null's 97.5th
  percentile. On the committed REF field the real gap is 249.361 survivors
  against a null of 0.198 (p97.5 0.203), 1258x. `verdict.json` records it as
  `noise_null`.
- G7b now compares the REF file's hash with the pin the committed
  `verdict.json` holds. The first run pinned; every later run compares.
- G8 now compares the heterogeneous arm's total with the uniform arm that
  `arm_pair` builds for the verdict at the primary level. Its control builds
  that arm from the median; the verdict would then move by 29%, and G8
  refuses it.
- A REFUSED-GATED run writes no `verdict.json`. Before this rule a gated run
  rewrote the G7b pin with the hash of the file it had just refused, and a
  second run of the same file passed. Found by the second review of
  2026-09-18.

Each gate has a planted defect that it must refuse, run by
`study2_response_layer/controls.py`: a flat noisy field (G6b), a different
real file against the pin (G7b), a uniform arm built from the median (G8).
All three are REFUSED-GATED, and `verdict.json` is byte-identical afterwards.

What does not change. The verdict, the four criteria and every number in
`verdict.json` before this amendment. The `2331x` figure the README once quoted
was `gap / sigma` printed under a 2-sigma label; the runner now prints
`gap / (2 sigma)`, which is the 1165x already stored as `ratio_to_2sigma`.
The noise-sigma mapping keeps its declared 1.645 divisor; the comment now
says what it is.
