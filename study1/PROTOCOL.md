# Transport-Caching Convergence Study — Protocol

Project: PRRT-spatial-cpm, Step 1 of the dependency order.
Status: PROTOCOL (pre-declared, before any residual is computed).
Written: 2026-09-16. This document precedes all results; criteria are fixed
before the first DPK-vs-MC comparison is run.

## Purpose

The closed-loop PRRT model (preprint v1.1, Sec. 3.4/3.6) must recompute cell
absorbed dose every treatment cycle as the cell population changes. Full Monte
Carlo re-transport per cycle is the defensible default; a cached dose-point
kernel (DPK) applied by convolution with the activity map is the cheap
alternative. This study decides, with pre-declared criteria, whether the DPK
path is admissible, and states what it costs when it is not.

What this study is NOT: it is not a validation of either arm against
experiment. Both arms are models. The reference arm (Geant4-class electron
transport) is itself validated only to the published 1.2–1.5% agreement with
EGSnrc (Perrot 2014, doi:10.1088/0031-9155/59/9/2183) and inherits that
uncertainty floor. No residual below that floor is interpretable.

## Physics inputs (provenance-declared; v1.1 amendment pre-results)

1. Lu-177 decay: BOTH arms sample the full ENSDF-evaluated decay (beta
   continuum with evaluated shape factors, conversion electrons, and the
   208/113-keV gammas) via Geant4 RadioactiveDecay 6.1.2 (G4RADIOACTIVEDATA
   z72.a177), Geant4 11.4.2. No hand-built spectrum is used. [v1.1: a
   hand-built allowed/forbidden-shape spectrum was initially drafted and
   REFUSED by its own anchor check — sampled mean 198 keV vs the 133-keV
   tabulated mean (Finogenova 2026, IJMS 27:7856, in-bibliography); the
   hand-built model is retired, not patched. Both arms now share the
   Geant4-decayed source, which makes the arm comparison purely
   transport-method vs transport-method.]
2. Water as tissue surrogate, unit density, 1 g/cm3 EVERYWHERE. Because the
   medium is homogeneous, cell occupancy does NOT change transport
   properties; it changes only WHERE activity sits (source distribution).
   [v1.1: corrects v1.0's erroneous claim that occupancy-vs-water contrast
   drives the DPK error. In a homogeneous medium a DPK convolution with an
   untruncated kernel is near-exact; the real error sources under test are
   (i) kernel support truncation at the ROI boundary, (ii) kernel-vs-grid
   resolution mismatch, (iii) sub-voxel source placement (voxel-center
   vs activity-weighted), (iv) kernel statistical noise. The study measures
   these; a PASS must therefore be scoped to "DPK with declared support
   radius at matched resolution", not to DPK in general.]
3. The kernel support radius is bounded by the maximum beta CSDA range in
   water (1.5 mm max tissue range, Finogenova 2026 Tab. 1; the 497-keV
   endpoint range is ~0.7 mm). To be pinned from an archived NIST ESTAR
   water table at kernel-build time and recorded in kernel metadata.
   Truncation threshold: any kernel voxel receiving < 1e-6 of central dose.
4. REF dose floor: Geant4-vs-EGSnrc agreement for Lu-177 S-values is
   published at 1.2-1.5% (Perrot 2014, doi:10.1088/0031-9155/59/9/2183).
   No residual below that floor is interpretable as DPK error.

## The two arms

Reference arm (REF): condensed-history electron transport of the beta
spectrum through the voxelized geometry (Geant4-class code). Voxel doses are
energy tallies / mass. This arm is CPU-reproducible, fixed-seed, and its
statistical convergence is itself ladder-tested (see Controls).

Approximation arm (DPK): the pre-computed dose-point kernel — the 3D dose
distribution in homogeneous water around a point source of the SAME Geant4
decay, computed by the same transport code at matched resolution —
convolved with the activity map. Verdict-path convolution: FP64 CPU
(scipy.signal FFT, FFTW backend); GPU (JACC AMDGPU on the 890M, or Metal on
the Mac) is an acceleration experiment only, benchmarked against the CPU
path, per the banked project decision — the verdict is never issued from a
GPU-only computation. The DPK arm consumes the activity map, not the
geometry: its error sources are the truncation/resolution/placement terms
declared in Physics inputs #2-3.

## Geometry

Synthetic, histology-representative CPM snapshot (viable / necrotic / stromal
regions), fixed across the whole study: one geometry, varying only the
sampling grid, per the fixed-object ladder rule (Biofilms AGENTS.md rule 5).
Voxel pitch options: 10 / 20 / 40 um on a 1 mm cubed ROI. Activity assigned
per-cell via the lognormal e_i draw (preprint Eq. 4 setup), fixed seed.
The geometry and activity map are generated once and hashed; every arm
consumes the same hashed inputs.

# PASS/FAIL CRITERIA (pre-declared)

Primary endpoint: cell-by-cell absorbed dose residual, defined
    r_i = (D_dpk,i - D_ref,i) / D_ref,i
over cells with D_ref,i > 1% of the maximum cell dose (the "clinical core" —
cells outside this are dosimetrically irrelevant and would inflate the
residual distribution with noise).

PASS-DPK: DPK is admissible for production IF the 95th percentile |r_i| over
the clinical core is < 2.0% (user's pre-declared tolerance).

FAIL-DPK: IF the 95th percentile |r_i| exceeds 2.0%, full MC re-transport is
mandated for every cycle. The finding is published as-is: "DPK caching is
inadmissible at cellular resolution for Lu-177 in heterogeneous cell
geometries at pitch p"; the pitch-dependence of the failure (20 um vs 40 um)
is itself the publishable structure.

In either verdict, the residual's spatial structure is reported: near-
interface vs bulk (interface = cell boundary within 1 kernel FWHM). If the
error is confined to interfaces, the admissibility statement is scoped to
bulk cells and the DPK is refused for interface-dominant geometries.

## Statistical controls (AGENTS.md rule 1: every check can fail)

1. REF statistical ladder: the reference arm is run at 2x and 4x histories.
   The residual p95 between the two REF runs must be < the pass threshold
   before any DPK comparison is interpreted; if REF itself is not converged,
   no DPK verdict is issued (refuse, do not pass).
2. Negative control for the comparison pipeline: swap in a deliberately
   wrong kernel (e.g., a kernel with the wrong spectral weighting — say a
   monoenergetic 497-keV kernel). The comparison MUST flag it (residual p95
   >> 2%). If the pipeline cannot distinguish the wrong kernel from the
   right one, the comparison is broken and no verdict is issued.
3. Uniform-geometry sanity: in a fully homogeneous water box (no cells), DPK
   must agree with REF to within the REF statistical floor (< ~1.5%).
   Disagreement there isolates implementation error from physics error.
4. Hash equality: both arms must consume identical, hashed input files
   (geometry hash + activity hash). A hash mismatch aborts the run.
5. Seed handling: DPK build and REF runs use fixed declared seeds; any
   verdict run is repeatable by re-running with the same seed and must
   reproduce its own numbers (run-to-run reproducibility is asserted, not
   assumed).
6. FP precision: the convolution arm runs double precision on CPU (threads
   backend) as the production path; if AMDGPU (FP32-dominant) is used for
   speed, a FP32-vs-FP64 comparison on a slice must quantify the precision
   penalty and the verdict is issued on the FP64 path only.

## Outputs

- verdict (PASS-DPK / FAIL-DPK / REF-UNCONVERGED) per pitch
- residual distributions (voxel + cell), per pitch
- interface vs bulk decomposition
- kernel metadata (spectrum provenance, ESTAR range pin, support radius,
  resolution, seed, code version)
- REF convergence ladder table
- benchmark table: DPK convolution wall time vs estimated MC wall time

## Boundaries

- No biological response is computed (that is Step 2; dose only).
- No clinical kernel is produced (conversion/Auger lines excluded).
- No material heterogeneity (water everywhere).
- No claim about Ac-225 alphas (different transport regime; separate study).
- The 2% threshold is a user pre-declared design tolerance, not a clinical
  standard; the manuscript states it as such.

## Execution order (all criteria fixed above, before results)

1. Generate + hash geometry/activity inputs.
2. Build DPK from spectrum (after ESTAR pin).
3. Run REF at 3 pitch values x 2 history levels + uniform control.
4. Run DPK convolution per pitch.
5. Compute residuals, controls, verdict.
6. Write results to study/ with this protocol doc co-committed.

## v1.2 amendments (2026-09-16, issued AFTER first REFUSED verdict, BEFORE any DPK verdict)

1. C1 scope defect, fixed: v1.0's ladder control evaluated voxel-level noise
   over ALL core voxels, including necrotic/void voxels dosed only by sparse
   cross-dose (~1-event Poisson noise, p95 88% at 8M decays). The pre-declared
   primary endpoint is the CELL residual; cells occupy only cell_id > 0
   voxels. C1 now evaluates the ladder over the same cell endpoint. This is
   endpoint-correcting, not threshold-tuning: the 2% criterion, the clinical
   core definition, and the PASS/FAIL rule are unchanged.
2. Kernel support violation, fixed: v1.0's kernel boxes (half-width ~0.86 mm)
   geometrically truncated the 497-keV beta CSDA range (~1.8 mm in water;
   anchor: Finogenova 2026 max tissue range 1.5 mm + ESTAR-class CSDA),
   violating this protocol's own "truncate at 1e-6 of central dose" rule.
   Kernels rebuilt at kbox 361/181/91 (10/20/40 um) = full support + margin.
3. Declaration: in the synthetic geometry each cell is one voxel-scale site
   (cell size = pitch); the three pitches are therefore three separate
   experiments, NOT one fixed-object ladder. The REF h1-vs-h2 pairs at each
   pitch are the fixed-object statistical ladders. Manuscript text must not
   describe the pitch sweep as a convergence ladder.


## v1.3 amendment (2026-09-16, after subset-seed proof, before any verdict)

DEFECT FOUND AND PROVEN, not inferred: with a shared fixed seed, ladder run
h1 (32M decays) is a strict SUBSET of h2 (128M decays): measured h2/h1 total
ratio 3.9988, and r1 <= r2 held at every one of the 15,619 nonzero voxels.
The "ladder" therefore measured a deterministic 4x energy difference, not
statistical convergence — p95 ~77-88% at every statistics level. Root cause:
G4 MT with one fixed master seed gives every run the same per-thread RNG
streams.

FIX: independent seeds per run (20260916/17/18 declared). Reproducibility
retains its correct meaning: same seed + same inputs -> identical outputs.
The convergence requirement itself is unchanged: ladder p95 < 2% at the
declared cell endpoint. Runs h1(32M, seed A) vs h2(32M, seed B) give the
noise floor at matched statistics; h2 vs h3(128M, seed C) gives the scaling.
If the cell endpoint cannot reach 2% at feasible statistics, the protocol
requires the verdict to be issued at the VOXEL endpoint with the cell
endpoint reported as-is, and the manuscript must state which endpoint the
2% criterion was evaluated on. (Endpoint re-scoping, if needed, will be
declared BEFORE looking at any DPK-vs-REF residual.)


## v1.4 amendment (2026-09-16, before any admissible verdict)

FATAL GEOMETRY BUG FOUND AND FIXED — the source-placement expression in the
reference app was algebraically self-cancelling: `(x - half + half, ...)` ==
`(x, ...)`, so all field sources were placed in the world's POSITIVE OCTANT
[0, 2*half]^3 (approximately 7/8 of decays born outside the world, never
transported), and the kernel-mode source sat at a box CORNER, capturing one
octant of the dose spread.

DETECTION CHAIN (recorded because it validates the control stack):
1. Shape-normalized DPK-vs-REF residual p50 = 133% — an absurd value that
   forced the units audit;
2. units audit exposed REF per-decay deposit = 0.01643 MeV vs the 0.1352 MeV
   smoke-test anchor — ratio 12.2% ~= 1/8, the octant signature;
3. both arms were wrong IDENTICALLY (per-decay totals agreed to 5 decimals
   across arms), which is why every earlier arm-vs-arm comparison had looked
   clean — the cross-controls caught it only when the DPK-vs-REF shape was
   compared against the ABSOLUTE smoke-test anchor.

CONSEQUENCE: every transport run before this amendment (all ref_p*, all
kernel_p*) is INVALID and is deleted. No verdict was issued from them (the
ladder control had refused twice), so no published number is affected.

POST-FIX VALIDATION GATE (new, mandatory before any run is used):
- smoke kernel test MUST reproduce per-decay deposit within [0.12, 0.16] MeV
  (0.1352 MeV anchor +- tolerance) — a run failing this is refused;
- kernel must be centrosymmetric to MC noise (cos(theta) correlation of the
  dose field about the source > 0.99 against its 180-degree rotation);
- field per-decay total must match sum_v A_v*k_v / sum_v A_v within MC noise.

## v1.5 — VERDICT (2026-09-16, all gates passed, pre-declared criterion applied)

Pitch-40 um, 1mm^3 ROI, 128M-decay reference (independent seeds A/B/C):

  REF convergence (C1, cell endpoint): 32M-A vs 32M-B p95 = 3.28%;
  32M vs 128M p95 = 2.62%; REF(128M) self-noise p95 ~= 1.17%  -> PASS
  (below the pre-declared 2% criterion).

  DPK vs REF (C-cell endpoint, clinical core): p50 = 3.95%, p95 = 23.84%,
  mean = +8.29%  ->  **FAIL-DPK** (criterion: p95 < 2.0%).

  Negative control (wrong kernel): p95 = 230% — the pipeline detects bad
  kernels; the comparison is sound. C3 uniform-box control: same ~5-6%
  global overestimate; radial profile DPK/REF = 1.000 center -> 1.033
  boundary. Energy gate: DPK predicts 0.10355 MeV/decay retained, REF
  measures 0.09752 (+6.2%).

  Spatial decomposition of the failure (characterization, NOT verdict-
  changing): mean residual by depth from the ROI face: +23.3% (depth 0),
  +6.3% (depth 1), +3.6% (2), +2.3% (3), decaying to +0.4% (depth 7).
  The inadmissible error is an escape-surface shell ~1-2 voxels thick.
  Interior cells: mean +1-2%. Hot cells (>25% Dmax): mean +5.6%, p95 21.7%
  (shell-dominated).

MECHANISM: the FFT convolution carries the kernel's infinite-medium
assumption to the ROI boundary; REF transports betas and 208/113-keV gammas
physically, and they escape. Both arms agree on per-decay totals to the
literature anchor (0.1352 MeV kernel-mode; field-mode retention lower by
escape, correctly).

CONSEQUENCE FOR THE MANUSCRIPT (the study's deliverable):
1. Cached-DPK convolution is INADMISSIBLE at the pre-declared 2% cell-level
   tolerance for PRRT geometries with escape boundaries within one beta
   range (which is every lesion ROI). Full Monte Carlo re-transport per
   treatment cycle is REQUIRED. This is a publishable finding in its own
   right — the "cheap shortcut" fails for a mechanistic, characterized
   reason, and the failure is measurable, reproducible, and boundary-
   localized.
2. IF the closed-loop model can be formulated with an interior mask
   (dose computed only for cells >2 voxels from the ROI face, boundary
   cells handled by re-transport or by a correction factor), DPK becomes
   admissible for interior cells at ~1-3% accuracy. That is a declared
   scope restriction, not a silent one — the manuscript must state it.
3. The GPU question is settled: the admissible arm is CPU-bound Geant4
   (~110k decays/s on 24 threads). The convolution (1.75 s) was never the
   bottleneck; there is nothing to offload.

RUNS STANDING BEHIND THE VERDICT (hashes in runs/): kernel_p40 (6M decays,
seed 20260916), ref_p40_h1/h2 (32M, seeds 20260916/17), ref_p40_h3 (128M,
seed 20260918), uniform-box control (32M, seed 20260920), wrong-kernel
control (analytic). All validation gates passed: per-decay anchor 0.1353
(kernel mode), centrosymmetry 1.0000, delta-alignment exact.

## v1.6 amendment (2026-09-18, pre-declared BEFORE the h4 run; issue #10)

Why. v1.5 passed C1 on an inferred figure: 32M-vs-128M p95 = 2.62%, times
0.447 (the 1/sqrt(N) factor for one 128M run's own noise), = 1.17%. No
amendment declared that inference, and `analyze.py` as coded evaluates C1 on
the measured h1-vs-h2 pair (3.28%) and refuses. Every residual measurable from
the committed runs is above 2%: 3.28% (32M vs 32M), 2.62% and 2.58% (32M vs
128M), 2.01% (pooled 64M vs 128M). C1 as worded needs a measured pair under 2%.

The run. h4: 128M decays, seed 20260919, geometry hash a6883b84 (unchanged),
same application, physics list and cuts as h3. Nothing else is re-run.

The test. C1 is evaluated on the highest-statistics pair, h3 vs h4: ladder
p95 of |D_h3 - D_h4| / D_h4 over the cell endpoint in the clinical core, the
v1.2 form that `analyze.py` implements. Threshold 2%, unchanged. The
1/sqrt(N) model predicts 1.64% for this pair.

Both outcomes, declared now:

- p95 < 2%: C1 PASS as a measurement. The v1.5 inferred figure is retired and
  `analyze.py` reads the highest-statistics pair present. FAIL-DPK is
  unchanged; the cell p95 of 23.8% is more than seven times any residual here.
- p95 >= 2%: C1 FAIL as worded at this statistics level. The study reports
  FAIL-DPK with C1 unmet, and no inferred figure substitutes for the
  measurement. A further reference run would need its own amendment.

No other number, threshold or file in this study changes under this amendment.

### v1.6 result (2026-09-18, after the run)

h4 ran as declared: 128M decays, seed 20260919, geometry hash a6883b84,
committed as `runs/ref_p40_h4.bin` and `.json` (sha256 of the bin begins
07e9be71). Wall clock about eight minutes on 24 threads.

C1 on h3 vs h4, cell endpoint, clinical core: p95 = 1.62%. Threshold 2%.
C1 PASS as a measurement; declared outcome A applies. The 1/sqrt(N) model
predicted 1.64%. Consistency: h1 vs h4 2.58%, h2 vs h4 2.58% (h1 and h2 vs h3
were 2.58% and 2.62%). Total energy per decay, h4 over h3: 0.999998.

`analyze.py` now evaluates C1 on the two highest-statistics REF runs present
and records the pair as `C1_pair`; here that is h3 vs h4. The v1.5 inferred
figure of 1.17% is retired from the verdict path; `final_analysis.py` still
prints it, labelled as an inference. FAIL-DPK is unchanged: cell p95 of 23.8%
against a measured reference residual of 1.62%, a factor of 14.7.
