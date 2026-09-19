# Closed-Loop Selection Study — Protocol

Project: PRRT-spatial-cpm, Step 3 of the dependency order.
Status: PROTOCOL (pre-declared, before any closed-loop run).
Written: 2026-09-17. This document precedes all results; the loop mechanics,
gates, and verdict criteria below were fixed before the first 4-cycle
simulation was executed. Builds on Study 1 (FAIL-DPK -> full re-transport
mandated) and Study 2 (PASS-DIVERGENCE -> response layer verified open-loop).

## Purpose

The closed-loop PRRT model (preprint v1.1, Sec. 3.6) treats each treatment
cycle as: transport on the CURRENT geometry/activity map -> protracted-LQ
response -> mitotic-catastrophe death at division entry -> repopulation ->
re-uptake proportional to the heritable SSTR2 expression trait. This study
executes that loop for a declared 4-cycle clinical regimen on the validated
Study-1/2 stack and verifies the selection dynamic: differential survival
under heterogeneous cross-dose must drive the population's expression
distribution {e_i} downward — the model-side reproduction of the
retreatment-failure observation (Feijtel 2021, Fransson 2024).

What this study is NOT: not a calibration of the CPM (no Hamiltonian
tuning; the full Potts surface-tension dynamics are Step 4+), not a
prediction for any patient lesion (synthetic geometry), and not a claim
about absolute TCP in clinical units.

## Fixed objects (hashed; consumed, never regenerated)

1. Initialization: study1_dpk_vs_mc/geometry_pitch40.npz,
   sha256 = a6883b84... (verified against geometry_pitch40_meta.json at
   every run start). T0 state: 13,978 cells (cell_id > 0) filling the
   viable band (11,000 sites) and the stromal band (2,978 sites); the 1,647
   necrotic voxels hold no cell. The tissue label has no computational role
   in this study: stromal cells share the expression draw, uptake and
   dynamics of viable cells (issue #14). Masks fixed for the whole study,
   lognormal expression
   (mu = 0, sigma = 0.6), activity map shape (CV = 0.658).
2. Physics engine: study1_dpk_vs_mc/build/study1_app (Geant4 11.4.2, MT,
   ENSDF RDM 6.1.2, the exact binary that passed the Study-1 validation
   gates: per-decay anchor 0.1353 MeV, centrosymmetry 1.0000, delta-
   alignment exact). Invocation contract:
   study1_app field 40.0 <geom_prefix> <out_prefix> <nEvents> [threads]
   [kbox] [seed]; consumes <geom_prefix>.bin (int32 cell_id[n^3] then
   float64 activity[n^3], C-order) + <geom_prefix>_meta.json; emits
   MeV/voxel C-order + sidecar JSON. DPK arm remains EXCLUDED.
3. Response engine: the locked Study-2 functions — SF_i =
   exp(-alpha*D_i - beta*G*D_i^2), alpha = 0.24 Gy^-1, beta = 0.06 Gy^-2,
   T_rep = 1.5 h, T = 96 h, G(96 h) = 0.04407 (Taylor-patched
   implementation, protocol v1.1 amendment 1), TCP endpoint not used
   here (the population IS the endpoint now; the stochastic realization
   replaces the Poisson expectation).

## Declared parameters (the clinical cycle; fixed before any run)

1. Regimen: 4 cycles, 8 weeks inter-cycle spacing (NETTER-1 convention),
   biological time; each cycle = 96 h active exposure (Step-2 declared
   T) at constant dose rate, then 8 weeks - 96 h of biological recovery
   and proliferation before the next cycle's transport.
2. Administered activity per cycle: A_admin fixed across all cycles.
   Scaling to clinical units is DECLARED, not fitted: the cycle-1 dose
   field is normalized so the mean viable-cell dose is D_mean = 10 Gy
   (Study-2 primary level); the same A_admin (in relative units) is
   re-applied every cycle. Per Choice B, total lesion activity per cycle
   is CONSTANT — if total dose dropped merely because cells died, TCP
   would fall for the wrong reason; the fixed-activity constraint
   isolates the spatial selection effect.
3. Uptake proportionality (preprint Eq. feedback, the load-bearing
   declared choice): per cycle, A_j = A_admin * e_i / sum_i(e_i) over
   viable cells i (one voxel per cell at pitch 40, Study-1 convention);
   non-viable voxels get zero activity. sum_j A_j = A_admin exactly.
4. Expression drift at division: e_daughter = e_parent * exp(eta),
   eta ~ N(0, sigma_div), sigma_div = 0.10 (declared standing variation;
   e stays strictly positive by construction).
5. Mitotic clock: doubling time T_d = 58.4 h (Tamborino 2025, NCI-H69);
   each surviving cell divides ONCE per inter-cycle window is NOT
   assumed — division attempts are Poisson with rate lambda_div =
   ln(2)/T_d per cell per biological day, realized per inter-cycle
   window; the population regrows toward the viable capacity.
6. Capacity: the T0 cellularized footprint (13,978 sites, viable plus
   stromal bands) is the carrying capacity. Division requires a free site;
   daughters inherit e with drift (declared #4). Necrotic voxels hold no
   cell and are never colonized. Stromal sites are occupied at T0 and
   refill after death like viable sites.
7. Death rule (preprint: death at mitosis entry): a cell that received
   cycle-k dose attempts its FIRST post-exposure division with survival
   probability SF_i; failure = mitotic catastrophe, site cleared. Cells
   surviving their first division are not re-dosed within the same
   cycle (declared simplification; intra-cycle re-dosing deferred).
8. CPM clock: 1 biological day = 100 Monte Carlo Steps (declared
   mapping; the Potts Hamiltonian is inactive in this study — the CPM
   provides the lattice bookkeeping, volume resorption of cleared sites,
   and the division/inheritance mechanics). 8 weeks = 56 days per
   inter-cycle window.
9. Transport statistics per cycle: 8M decays (justified: Study-1 ladder
   showed cell-endpoint noise ~ sqrt(128/8) x 1.17% ~ 4.7% p95 at 8M —
   larger than the per-cycle selection signal per cell, but the endpoint
   is the POPULATION mean over 13,978 cells, where the noise averages
   down; the verdict criterion is population-level, declared in
   Endpoints). Independent transport seed per cycle per arm (Study-1
   v1.3 lesson: shared seeds make runs non-independent).
10. Population-dynamics seed: separate seed for the death/division/
    drift realization; closed-loop and open-loop arms each get 5
    independent dynamics seeds (declared below).

## Iteration sequence (one cycle; per Choice B)

1. TRANSPORT: write the current activity map to <geom_prefix>.bin,
   invoke study1_app (field mode, 8M decays, declared seed), read the
   MeV/voxel field.
2. RESPONSE: cell doses D_i from the field (Study-2 declared scaling,
   D_mean = 10 Gy at cycle-1 shape; same physical conversion every
   cycle); SF_i from the locked LQ functions.
3. GROWTH/DEATH: realize first-division survival per cell (Bernoulli,
   SF_i); catastrophied cells cleared, sites freed; survivors divide
   per the Poisson mitotic clock into free sites with expression
   inheritance (#4) until capacity or window end.
4. RE-UPTAKE: distribute the fixed A_admin over the new viable map
   proportional to e_i (declared #3). Cycle counter ++.

## The population-dynamics resolution (declared, pre-run)

Deaths are NOT cumulative across cycles without repopulation: with the
declared doubling time, the ~1,220 cycle-1 survivors (13,978 x mean SF
= 0.0875 at 10 Gy) repopulate the 13,978-site cellularized footprint in ~3.5
doublings (~8.5 days), far inside the 8-week window. Every cycle
therefore STARTS at capacity with ~10 Gy mean dose; the sawtooth
(plunge, rebound) is per-cycle, and the selection signal is carried by
the survivor pool's trait distribution. This resolution is what makes
the fixed-activity constraint self-consistent across 4 cycles.

## Negative controls and gates (every check can fail)

K1 — Activity conservation: |sum_j A_j - A_admin| / A_admin <= 1e-12 at
     EVERY cycle start, including cycles 2-4 AFTER deaths and
     re-uptake. Enforced at write time; violation aborts the run.
K2 — Null selection (no-selection gate): uniform {e_i}, sigma_div = 0,
     uniform dose -> population mean e must be bit-stable across all 4
     cycles (all SF_i identical -> all death/division decisions
     trait-independent; any drift exposes trait-biased bookkeeping).
K3 — Extinction gate: dose scale x1000 (1000 Gy/cycle) must clear the
     lattice at cycle 1; the driver must handle an empty activity map
     gracefully (zero-activity transport run returns a zero dose field
     without error; population stays 0; run terminates with EXTINCT
     status, not a crash or NaN).
K4 — Open-loop control (Choice A): 4 cycles with the activity map
     FROZEN at cycle 1's distribution (re-applied blindly, same total
     activity, ignoring deaths). PRE-DECLARED PREDICTION: open-loop
     mean e does not drift (death is trait-independent in this arm by
     construction — the dose each cell sees is fixed, so SF_i is
     fixed, and daughters inherit e with symmetric drift); any open-
     loop drift > 0 (in either direction, beyond seed spread) falsifies
     the driver. The closed-vs-open CONTRAST is the selection signal.

## PASS/FAIL CRITERIA (pre-declared)

Endpoints (computed on ln e; the lognormal trait makes the log the
natural scale; geometric-mean drift is the reported quantity):

E1 (primary): closed-loop mean ln e at cycle-4 start < closed-loop mean
   ln e at cycle-1 start, AND the closed-vs-open contrast
   Delta_sel = (mean ln e closed, cyc4) - (mean ln e open, cyc4) is
   negative beyond the open-loop seed spread: Delta_closed_mean <
   min over the 5 open-loop replicates. Both arms run 5 independent
   dynamics seeds; the closed-loop seed spread is reported as-is.
E2 trajectory: viable cell count shows the declared per-cycle sawtooth;
   AND the per-cycle survivor count RISES across cycles under fixed
   A_admin (cycle-2+ survivors > cycle-1 survivors) — the pre-declared
   quantitative signature of selection (a more resistant population
   loses fewer cells to the same dose). A monotonically falling
   survivor count across cycles 1-3 falsifies the selection claim even
   if E1 passes spuriously.
E3 mechanism check (reported, not gated): per-cycle dose CV over
   viable cells; the manuscript states whether selection saturates via
   cross-dose flattening.

PASS-SELECTION: E1 (both parts) AND E2 (both parts).
FAIL-SELECTION: any part fails; reported with the failing part.
REFUSED-GATED: any gate K1-K4 fails; no verdict issued.

## Outputs

- verdict + per-gate results with measured values;
- per-cycle: viable count, mean/p50/p95 of e (and ln e), dose mean/CV
  over viable cells, activity CV, survivor counts;
- expression distribution snapshots (cycle 1 vs cycle 4, both arms);
- closed-vs-open contrast with the 5-seed spreads;
- all transport run metadata (seeds, decays, app sha).

## Boundaries

- Full Potts Hamiltonian (lambda_V, lambda_S, J) inactive: cells occupy
  fixed sites; only bookkeeping-level CPM mechanics (clearing, filling,
  inheritance) run. Declared: Step 4.
- No oxygenation, no phase-dependent radiosensitivity, no immune
  clearance, no stroma dynamics.
- Single lesion, one geometry, one regimen; scans are Step 4.
- Absolute clinical dose calibration is declared (D_mean = 10 Gy), not
  derived from a pharmacokinetic model; the loop's STRUCTURE (selection
  direction) is the claim, magnitudes are conditional on the declared
  regimen.
- 8M-decay transport per cycle: per-cell dose noise ~4.7% p95 is
  accepted and its propagation through SF is covered by the 5-seed
  dynamics spread; a 128M-decay confirmation of the final cycle-4
  field is declared as the robustness run IF the verdict is within a
  factor 2 of any criterion boundary.

## Execution order (all criteria fixed above, before results)

1. Verify fixed-object hashes; smoke-test the app invocation (1M decays
   on the T0 map reproduces the Study-1 field's per-viable-cell MEAN
   shape to within the declared statistical tolerance at 8M vs 128M
   scaling — coarse check, 10% on the mean, since the cycle-1 dose
   field must reproduce the Study-2 foundation).
2. Gates K1-K4 (K2/K3 on short 1-cycle runs; K4 on the full 4-cycle
   schedule).
3. Full runs: closed-loop x5 dynamics seeds, open-loop x5, extinction.
4. Endpoints E1-E3; verdict.
5. Write verdict.json + this protocol co-committed with the driver.

## v1.1 amendment (2026-09-17, issued during driver design, BEFORE any
## closed-loop or open-loop run; no endpoint has been computed)

K4's pre-declared prediction was WRONG as written, caught in pre-run
design review: v1.0 predicted "open-loop mean e does not drift; any
open-loop drift falsifies the driver." That prediction ignores the
arm's own construction: the frozen activity map is cycle-1's map,
which is e-CORRELATED (activity ∝ e_i). Re-applying it blindly keeps
killing the occupants of high-e sites in every cycle, so the open-loop
arm ALSO drifts downward — it just cannot COMPOUND the drift, because
re-uptake never re-concentrates activity on the enriched survivor
pool. The pre-declared prediction is corrected to:

  0 <= |Delta_open| < |Delta_closed| — the open-loop arm drifts
  weakly (frozen-site selection only), the closed-loop arm drifts
  strictly more (compounding selection), and the CONTRAST
  Delta_sel = |Delta_closed_mean| - |Delta_open_mean| > 0 with
  Delta_closed_mean below the open-loop seed spread in the direction
  of selection is the selection signal. An open-loop drift >= the
  closed-loop drift falsifies the driver. The K4 falsifier "any
  open-loop drift" is retired; it was wrong on the arm's own
  mechanics. E1's contrast definition is updated accordingly:
  Delta_sel = |Delta_closed| - |Delta_open| (population-mean ln e,
  cycle 4 vs cycle 1) must be positive, and Delta_closed_mean must
  exceed the worst (least-drifted) open-loop replicate.

No other criterion, gate, or parameter is changed. This amendment is
endpoint-correcting (the direction of the expected open-loop drift),
issued before either arm produced a number.

## v1.2 amendment (2026-09-17, issued during driver design, BEFORE any
## closed-loop or open-loop run; no endpoint computed)

E2's second clause was WRONG as written, caught by deriving the Jensen
arithmetic under this protocol's own declared dose constraint. v1.0
gated on "per-cycle survivor count RISES across cycles (the selection
signature)". Under the fixed-total-activity constraint (declared #2)
the sign is the opposite:

  Selection enriches survivors in low-e; re-uptake renormalizes
  activity ∝ e over the survivors, so the cycle-2 activity pattern is
  the T0 pattern RESTRICTED to the low tail and re-normalized — a
  COMPRESSED pattern (ln-e sigma among refilled cells ~ 0.19-0.25 vs
  T0 0.6, from ~3.5 lineage doublings of sigma_div = 0.1 drift). A
  compressed activity pattern flattens the dose distribution; by the
  same Jensen convexity that Step 2 established, E[SF] - SF(mean)
  shrinks with the variance, so per-cycle survivors FALL from the
  cycle-1 value (~0.0875 x CAP) toward the uniform-dose value
  (~0.0696 x CAP) and plateau. Falling per-cycle survivors under a
  fixed administered activity are the EXPECTED consequence of
  selection flattening the dose map — less dose wasted on hot cells,
  not a weakening of the treatment. The retreatment-failure signal
  this study tests is the RATCHET of the expression distribution
  (E1) plus population persistence (sawtooth refill to capacity),
  NOT a rising per-cycle survivor count.

Amendment: E2 is re-scoped to (i) sawtooth trajectory — deaths at
every cycle, refill to capacity before the next cycle (n_start = CAP
at every cycle start in the closed arm), population never extinguished
within the 4 cycles; and (ii) survivor-count direction ACROSS CYCLES
is RECLASSIFIED from gate to REPORTED DIAGNOSTIC, with the two
competing effects named (composition selection raises mean SF;
variance compression lowers the Jensen excess; net sign is an
empirical output and is reported with the dose-CV trajectory, E3).

Open-loop refinement (mechanics, not a new gate): with frozen
slot-weights and well-mixed refill, occupancy at cycle-2 start is
e-independent of site, so the frozen map's deaths are e-independent
from cycle 2 on — the open-loop arm's drift is a ONE-SHOT shift
(cycle 1's e-correlated selection) that plateaus. The closed-loop arm
compounds every cycle. Declared prediction, sharper than v1.1:
Delta_open ~= the single-cycle selection differential (all replicates
<= 0, similar magnitude), Delta_closed strictly more negative and
growing per cycle. E1's criterion (closed below the worst open
replicate) unchanged.

Transport-economy declaration: the open arm's cycles 2-4 use the
frozen T0 activity map with a fixed transport seed — bit-identical
inputs — so ONE shared cycle-1 transport run serves all 10 replicates
and the extinction arm (x100 dose scale on the same field). The
closed arm transports per (dynamics seed, cycle >= 2): 15 runs.
Total: 16 transport runs at 8M decays.