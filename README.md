# PRRT-spatial-CPM: closing the loop from measured activity to cell-resolved response

Code, protocols, and verdicts for three executed studies behind the preprint
*Closing the Loop from Measured Activity to Cell-Resolved Response*
(`preprint/`). The framework couples synthetic, histology-representative
Cellular Potts geometry to Lu-177 radiotherapy transport and a heritable SSTR2
expression trait, and tests whether retreatment failure in neuroendocrine
tumors is an emergent property of spatial heterogeneity. The executed studies
use a synthetic concentric-shell lesion, and the Potts Hamiltonian is inactive
in them (`study3_closed_loop/PROTOCOL.md`).

The framework extends the one-way transport foundation of
[Kinder & Faulkner 2026, Biofilms](https://github.com/ffinkdevs/Biofilms)
(the `coupling/` HDF5 snapshot -> transport -> `import_dose_field!`
architecture) to the neuroendocrine, CPM, histological realm.

## The three studies (each pre-declared, gated, committed with its verdict)

| Study | Question | Verdict |
|---|---|---|
| 1 — transport caching | Is a cached dose-point kernel (FFT convolution) admissible vs full Geant4 re-transport? | **FAIL-DPK**: cell p95 = 23.8% (criterion 2%); boundary-shell mechanism; full re-transport mandated |
| 2 — open-loop response | Do protracted LQ kinetics on the validated dose field reproduce the Mellhammar TCP divergence? | **PASS-DIVERGENCE**: TCP_het < TCP_unif at every dose level; Jensen mechanism; 1165x the 2-sigma noise bar |
| 3 — closed loop | Does re-uptake ∝ heritable expression drive clone selection over 4 clinical cycles? | **PASS-SELECTION**: ln-e drift −0.484 closed vs −0.189 open-loop control; contrast +0.295; compounding only when the loop closes |

Every protocol was written and amended **before** any verdict-issuing run,
with negative-control gates that can fail (and did: the gate ladder refused
three defective G(T) implementations in Study 2 and two of the protocol's
own pre-declared predictions were corrected by amendment before the runs —
see each `PROTOCOL.md` amendment history).

## Repository layout

```
preprint/                 manuscript (LaTeX + PDF + bibliography data)
study1/                   DPK vs Monte Carlo transport study, Geant4 app
                          source (src/, CMakeLists.txt) and pitch-40 runs
study2_response_layer/    open-loop LQ + Lea-Catcheside TCP verification
study3_closed_loop/       4-cycle closed-loop selection with re-transport
figures/                  Study-3 summary figures
paraview/                 ParaView volumes: cycle-1 vs cycle-4 .vti + .pvd
                          series + the exporter that made them
data/                     hashed T0 geometry + cycle-4 state export
spec/                     response.growth_survival module spec + test vectors
```

## Implementing the response layer elsewhere

`spec/response_growth_survival.md` is the specification for the
`response.growth_survival` dividing/dying module, written for ports that add
the layer to their own CPM. It pins the constants, gives the test vectors, and
declares the two acceptable death/refill semantics.

```bash
python3 spec/check_vectors.py --controls     # vectors + gates G-N, G-O, G-P, G-C
python3 spec/check_selection.py --controls   # six population gates, reduced lattice
```

The first asserts every vector in the spec and runs the four closed-form gates.
The second runs six population gates on 2000 sites and 5 paired seeds. G-S and
G-D come from the spec. G-Q, G-M, G-B and G-H come from its amendments.

`--controls` feeds each gate an implementation it must reject, so a gate that
has stopped checking anything is visible rather than silently green. Both
scripts use the standard library only. The second needs Python 3.12 or later
for `random.binomialvariate`.

`spec/pr_checklist.md` is the acceptance criteria for a port's pull request. It
states the declarations the PR description must make and the tests it must
ship. It also names four hazards read from the Odin port at
`ffinkdevs/Biofilms@4788dcd1`, each with its status from the review of `f3bf6d6`.

## Reproducing

Requirements: Python 3.12+ (numpy, scipy), Geant4 11.4.2 with
RadioactiveDecay 6.1.2 (ENSDF), 24-thread CPU recommended.

```bash
# Transport app, used by Studies 1 and 3. It derives from the Geant4 example
# examples/extended/radioactivedecay/rdecay02, which is not copied here.
source study3_closed_loop/g4env.sh      # export GEANT4_PREFIX first if yours differs
cmake -S study1 -B study1/build -DCMAKE_BUILD_TYPE=Release
cmake --build study1/build -j           # produces study1/build/study1_app

# Study 1 analysis (seconds; reads the committed pitch-40 runs in study1/runs/)
python study1/final_analysis.py

# Study 2 (minutes; needs the 128M-decay reference run, study1/runs/ref_p40_h3.bin)
python study2_response_layer/study2_run.py

# Study 3 (16 transport runs x ~75 s + dynamics). Every mode of study3_run.py,
# the default `gates` included, invokes study1/build/study1_app.
python study3_closed_loop/study3_run.py full
```

`study1/runs/` holds the pitch-40 reference, kernel and DPK outputs the
analysis scripts read. `study3_closed_loop/runs/` holds the four per-cycle dose
fields and the cycle-4 activity map of the representative closed-loop
replicate, which the figure and ParaView exporters read. Pitch-10 and pitch-20
runs were never made, and `study1/analyze.py` reports them as missing input.

The hashed fixed objects (`data/geometry_pitch40.npz`,
sha256 `a6883b84...`) pin the synthetic lesion geometry: 25^3 voxels at
40 um pitch, 13,978 cells across the viable and stromal bands (11,000 and
2,978; the tissue label has no computational role), lognormal expression
(sigma = 0.6), declared masks (viable/necrotic/stroma).

## Declared parameters (not fitted)

Radiobiology from [Tamborino et al. 2025](https://doi.org/10.2967/jnumed.125.269470)
(EBRT intrinsic fit for NCI-H69): alpha = 0.24 Gy^-1, beta = 0.06 Gy^-2;
repair half-time T_rep = 1.5 h (declared prior, consistency-gated against
their measured 1/G ~ 12); 96 h active exposure per cycle; doubling time
58.4 h; expression drift sigma_div = 0.10 per division; 4 cycles x 8 weeks
(NETTER-1 convention); mean viable dose 10 Gy/cycle (declared scaling).

## ParaView

`paraview/prrt_c1.vti` and `paraview/prrt_c4.vti` are VTK ImageData volumes
of the representative closed-loop run (cycle 1 and cycle 4), with point
arrays `cell_id`, `tissue_type`, `activity_relative` (K1-normalized uptake),
`dose_Gy`, `SF`, `ln_expression`, and a `.pvd` series
(`paraview/prrt_cycles.pvd`). Spacing is the declared 40 um per site with
its source recorded in the file's field data, following the evidence
discipline of the Biofilms `export_vti.jl` exporter. Axis convention is
pinned by an embedded orientation probe (x-fastest points, X=i, Y=j, Z=k).
Regenerate with `pvpython paraview/export_vti_prrt.py` (ParaView 6.2 /
VTK 9.7).

`paraview/hist4d_c1.vti` to `hist4d_c4.vti` and `paraview/histology_4d.pvd` are
four-cycle expression-state volumes on the same synthetic geometry. No
histological image is an input to them.

## License

CC0-1.0 (matching the Biofilms upstream). The preprint text is the
citation of record for the specification.