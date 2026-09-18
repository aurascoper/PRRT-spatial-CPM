# ISSUE SPEC — `response.growth_survival`: the dividing/dying CPM module

**Status: READY FOR IMPLEMENTATION.** The complete specification follows:
state variables, pinned constants, test vectors, transition rules, semantics
choices, and gates.

Reference implementations are linked at the bottom. Everything here is derived
from two executed, committed, gate-verified studies. Nothing is invented for
this issue.

> **Provenance.** The text below was written for
> [aurascoper/PRRT-spatial-CPM#1](https://github.com/aurascoper/PRRT-spatial-CPM/issues/1)
> and circulated as `ffinkdevs_growth_survival_SPEC.md`. It is committed here
> unchanged so that section numbers stay citable. Corrections found after
> circulation are listed at the end under **Corrections**, and are not folded
> silently into the body.
>
> Two checkers accompany it. `spec/check_vectors.py` asserts every vector in
> section 3 and runs the closed-form gates G-N, G-O, G-P and G-C.
> `spec/check_selection.py` runs the population gates G-S and G-D on a reduced
> lattice. Both take `--controls`, which feeds every gate an input it must
> reject.

## 0. The one-paragraph summary

A cell that has accumulated dose D during the exposure window attempts division
at its next mitosis.

It survives with probability `SF = exp(-alpha*D - beta*G*D^2)`. If it survives,
it divides normally and passes its traits (with small heritable drift) to
daughters. If it fails, it dies by mitotic catastrophe and its lattice site is
cleared and refilled. G is the Lea-Catcheside protraction factor, already
computed for the exposure window; it is a constant per cycle, not a per-cell
computation.

*(Correction 1 applies to the refill mechanism named in this paragraph.)*

## 1. State variables to add

| Variable | Type | Update rule |
|---|---|---|
| `accumulated_dose` | Float64, per cell | Set once per cycle from the transport/dose-accumulation step. Not decremented. |
| `cell_state` | Int8 enum | 1 = viable, 2 = mitotic (transient within the division step), 3 = dying/necrotic |
| `expression_e` | Float64, per cell | Heritable trait; daughters inherit parent value `* exp(sigma_div * N(0,1))` |

Existing volume/target-volume machinery is reused unchanged. No new Hamiltonian
terms are required for the minimal module.

`sigma_div = 0.10`. See **Answer 1** for its provenance; it was absent from this
section when the spec circulated.

## 2. Pinned constants (do not tune; declared provenance)

- `alpha = 0.24` Gy^-1 (0.24 ± 0.05)
- `beta = 0.06` Gy^-2 (0.06 ± 0.02)
- Source: EBRT fit on the NCI-H69 neuroendocrine line, Tamborino et al. 2025
  (JNM 66:1291). The PRRT-fitted coefficients from the same table are NOT used:
  they already embed protraction, and applying G on top of them double-counts
  it.
- `T_rep = 1.5` h (repair half-time) → `mu = ln(2)/T_rep = 0.462098` /h.
  Declared modelling choice; Tamborino publishes no explicit `T_rep`.
- `G(T)` for the exposure window T. Write `x = mu*T`. The exact form is
  `G = 2*(x - 1 + exp(-x))/x^2`. The true small-x limit is `G -> 1`: maximum
  interaction. The numerator `x - 1 + exp(-x) = x^2/2 - x^3/6 + ...` collapses
  to 0/0 in float64.
  - if `x < 1e-3`: `G = 1 - x/3 + x^2/12` (Taylor branch; the limit is 1)
  - WARNING: a form sometimes transcribed, `0.5*x*(1 - x/3)`, is the Taylor
    expansion of a different quantity and dives to 0. Reject it.
  - else: `G = 2*(x + expm1(-x))/x^2` (standard branch)
  - Use `expm1`, never a hand-written `exp(-x)`: the identity
    `x - 1 + exp(-x) == x + expm1(-x)` is exact, and `expm1` keeps full
    relative accuracy.
  - The naive form returns `G(1e-9 h) = 0.0` instead of 1.0. That silent
    failure is what gate G5c caught in the reference study. The G-P gate in
    section 5 catches it in your port.

## 3. Pre-computed test vectors (assert these exactly)

G values: `G(24 h) = 0.164076`, `G(96 h) = 0.044068`. Small-x branch vectors:
`G(T = 1e-9 h) = 0.999999999846`, `G(T = 1e-6 h) = 0.999999845967`. The branch
is continuous at `x = 1e-3`. *(Correction 2 gives the tolerance that applies
there.)*

SF at `G(96 h)`, alpha/beta pinned above:

| D (Gy) | SF (float64) |
|---|---|
| 0 | 1.000000 |
| 2 | 6.122734e-01 |
| 5 | 2.819285e-01 |
| 10 | 6.964060e-02 |
| 20 | 2.858008e-03 |
| 40 | 9.851017e-07 |
| 1000 | 0.0 (underflow; must not NaN) |

Protraction direction check at D = 10 Gy: `SF(T=24 h) = 0.033896 <
SF(T=96 h) = 0.069641`. Shorter exposure at fixed dose is strictly more lethal
(higher dose rate, less repair). If your implementation reverses this, the G
evaluator is broken.

Tolerance: assert to 6 significant digits; both reference implementations agree
to that level.

Machine-readable copy: `spec/reference_vectors.json`.

## 4. The stochastic transition (the coin flip)

At the first post-exposure division attempt of each cell:

1. Draw `RNG` uniform on [0, 1).
2. If `RNG <= SF_i`: the cell divides. Spawn the daughter per the existing
   division mechanics; halve volumes per the CPM convention; daughters inherit
   `expression_e` with drift `e_parent * exp(sigma_div * N(0,1))`.
3. If `RNG > SF_i`: set `cell_state = dying`. Death is by mitotic catastrophe;
   the cell does not divide.

Declared semantics (the PR must declare which was implemented; both are
gate-testable):

- **Semantics A (reference implementation):** clear the site immediately after
  the failed draw. Refill to capacity during the same Monte Carlo window, with
  heritable trait drift. *(Correction 1: the placement rule stated here was
  wrong.)*
- **Semantics B (CPM-native alternative):** set `V_target = 0` for the dying
  cell and let Metropolis dynamics resorb it over subsequent steps. Refill then
  happens through normal neighbour competition.

Either is acceptable. Mixing them silently is not.

## 5. Gates (the PR fails review if any of these cannot fail)

- **G-N (null gate):** D = 0 ⇒ SF = 1.0 exactly. A cell receiving zero dose
  must never die by this module.
- **G-O (obliteration gate):** D = 1000 Gy ⇒ SF underflows to 0.0 in float64.
  Must be exactly 0.0, never NaN, never negative.
- **G-P (protraction gate):** at fixed D, `SF(T=24 h) < SF(T=96 h)`, strictly,
  at every tested dose level. Also: G at tiny exposures must approach 1, not 0
  (see the small-x vectors in section 3). A G evaluator that returns 0.0 at
  `x = 1e-9` has the naive-cancellation defect. The Taylor branch in section 2
  is the only accepted fix.
- **G-C (consistency gate):** `1/G(96 h)` with `T_rep = 1.5 h` must fall in the
  range the Tamborino PRRT fits imply (about 12 to 23). The declared value is
  22.7. A wild miss means `T_rep` or the G formula is wrong.
- **G-S (selection-null gate):** uniform expression, zero drift, zero dose ⇒
  zero net drift in mean expression across any number of cycles. If mean e
  shifts with no dose gradient, the death/refill mechanics are biasing the
  trait. That is a bug, not a finding.
- **G-D (direction gate):** under heterogeneous dose, mean expression must FALL
  over cycles (selection ratchet). In the verified implementation the
  closed-loop fall is -0.48 log-units over 4 cycles against -0.19 for the
  frozen-map open-loop control. Your implementation need not reproduce those
  exact numbers (different geometry), but the sign and the closed > open
  ordering are the physics.

### Where these gates live in this repository

The spec names gates `G-N` through `G-D`. The committed studies predate that
naming and use their own. The map:

| Spec gate | Committed as | File | Assertion as coded |
|---|---|---|---|
| G-N | `G1 null` | `study2_response_layer/study2_run.py:118` | `SF(0 Gy) == 1.0`, bit-exact |
| G-O | `G2 obliteration` | `study2_response_layer/study2_run.py:120` | `SF(1000 Gy) <= 1e-100` |
| G-P | `G3 rate/protraction` | `study2_response_layer/study2_run.py:129` | `SF(D,24h) < SF(D,96h)` at D ∈ {2,8,10} |
| G-P (small-x half) | `G5c G(T->0) -> 1` | `study2_response_layer/study2_run.py:109` | `abs(G(0) - 1.0) < 1e-9` |
| G-C | `G5d Tamborino consistency` | `study2_response_layer/study2_run.py:112` | `1.0 < 1/G(96h) < 100.0` |
| G-S | `K2 null-selection` | `study3_closed_loop/study3_run.py:280` | mean e `== 1.0` exactly, 4 cycles |
| G-D | endpoints `E1a`/`E1b` | `study3_closed_loop/study3_run.py:341-347` | closed drift below worst open replicate; contrast > 0 |

`spec/check_selection.py` runs G-S and G-D on a reduced lattice, 2000 sites and
5 paired seeds, so both are executable without a transport code. The absolute
drift there is not the study's and is not meant to be.

Three further gates come from the amendments below and have no counterpart in
the committed studies: `G-Q` (arrested cells still die), `G-B` (dose is not
banked across cycles) and `G-H` (no cell hoards the activity).

`spec/pr_checklist.md` states what a port's pull request must declare, which
tests it must ship, and four hazards specific to the Odin port.

Two differences worth knowing before you write assertions:

- **G-O is stricter here than in the committed gate.** The spec demands exactly
  `0.0`; `G2` demands `<= 1e-100`. Both hold: `SF(1000 Gy)` underflows to
  exactly 0.0. Assert the spec's version.
- **G-C is narrower here than in the committed gate.** The spec band is 12 to
  23; `G5d` allows 1 to 100. The declared value 22.6922 is inside both, but it
  is 1.3% below the spec's upper edge. A port that changes `T_rep` even
  slightly upward will breach the spec band while still passing `G5d`.
- **G-D is not a `gate()` call in study 3.** It is a verdict criterion. There is
  no `gate("K4", ...)` at all; `K4` is declared in `PROTOCOL.md` and its
  original falsifier was retired by the v1.1 amendment.

## 6. Reference implementations (the math, working and verified)

- `study2_response_layer/study2_run.py` — protracted LQ survival with the
  Taylor-patched G evaluator; the G-N/G-O/G-P/G-C equivalents are enforced
  there under the names in the table above.
- `study3_closed_loop/study3_run.py` — the mitotic-catastrophe death plus
  refill loop, per-(cycle, replicate) seeding. *(Correction 1: this loop is not
  neighbour expansion.)*

Both are committed with protocols and verdicts in `aurascoper/PRRT-spatial-CPM`.

## 7. Scope boundary (what this module is NOT)

No bystander signalling. No immune response. No continuous-dose-rate
integration during the cycle (dose arrives as a per-cycle accumulated field).
No clonogen expansion beyond the existing refill mechanics.

The module is the minimal dividing/dying extension that flips
`response.growth_survival` from `unsupported_by_current_model` to supported,
nothing more. Every addition beyond this spec needs its own declared provenance
in the PR description.

---

# Corrections

Found after the spec circulated. Listed, not folded into the body, so that a
reader who implemented against the circulated text can see what moved.

## Correction 1 — the reference refill is well-mixed, not neighbour expansion

**2026-09-17.** Sections 0, 4 and 6 described the Semantics A refill as
"neighbour expansion". The reference implementation places daughters uniformly
at random over every free site in the lesion, with no spatial locality at all.

`study3_closed_loop/study3_run.py:163`:

```python
tgt = free_idx[rng.permutation(free_idx.size)[:k]]
```

`free_idx` is every unoccupied slot in the whole viable footprint, not the
dying cell's neighbourhood. The file says so itself at line 14: `well-mixed
placement; Potts neighborhood dynamics are Step 4, per protocol Boundaries`.
The string `neighbor` appears nowhere else in that file.

**What this changes.** Semantics A is *well-mixed placement*, and it is a
declared stand-in, because study 3 ran without a Potts Hamiltonian. Semantics B
is the one the manuscript specifies. A port choosing A to match the reference
would reproduce the stand-in, not the model.

## Correction 2 — the branch-continuity tolerance is 2e-11, not 1e-11

**2026-09-17.** Section 3 stated that the Taylor and standard branches "agree
to 1e-11" at `x = 1e-3`. They agree to 1.686e-11.

The gap is the first dropped Taylor term, `x^3/60 = 1.67e-11` at `x = 1e-3`. It
is mathematics, not float noise, so a literal test at 1e-11 fails on a correct
implementation. `spec/reference_vectors.json` records a tolerance of 5e-11.

## Correction 3 — G-S as worded cannot fail

**2026-09-17.** Section 5 states G-S as "uniform expression, zero drift, zero
dose". At zero dose every cell has `SF = 1.0`, so no cell dies, so no refill
runs. The death and refill mechanics the gate exists to check never execute.

The committed arm uses uniform dose, not zero dose.
`study3_closed_loop/study3_run.py:253-256`:

```python
def gate_k2():
    """Null selection: uniform e, sigma_div=0, uniform dose, 4 cycles."""
```

Uniform 10 Gy kills most of the population and the survivors refill, so the
mechanics do run and a biased refill is visible. `spec/check_selection.py` runs
both forms and prints the deaths per cycle for each, so the vacuous one is
visible rather than silently green:

```
zero dose:    deaths/cycle [0, 0, 0, 0]
uniform dose: deaths/cycle [1870, 1895, 1875, 1866]
```

Implement the uniform-dose form. Keep the zero-dose form too if you like; it is
cheap, and it does catch a `SF` that is wrong at `D = 0`.

*(Correction 6 narrows what the uniform-dose form can detect.)*

## Correction 4 — the checker enforced five significant digits while claiming six

**2026-09-18. Found by the Odin port**, which bit-matched Julia on the vectors
and reported a rounding problem.

`spec/check_vectors.py` computed its tolerance as `10**(-sig + 1)`. At `sig = 6`
the formula gives 1e-5, a full decade looser than the 5e-6 that six significant
digits means.

The vectors file hid the gap. `SF_at_G96` was stored in `%.6e`, which gives seven
significant digits. `G_of_T` and `protraction_at_10Gy` were stored in `%.6f`,
which gives five for any value below 1. One `_tolerance` string claimed six for
all of them.

`SF(10 Gy, 96 h)` was stored as 0.069641, against a true 0.0696405992270619. The
relative error is 5.75e-6, so it failed six digits and passed only on the slack.
The checker also printed `round(got, 6)` on that line, so the failure displayed
as two identical numbers.

**What changed.** The tolerance is now `5 * 10**(-sig)`. Every float vector is
stored at seven significant digits and asserted at six, leaving headroom. The
display prints full precision.

A new control proves the tolerance refuses the old five-digit value, accepts the
seven-digit one, and that the old formula accepted the five-digit value. Without
the third case nothing would show the change did anything.

The file also claimed `_generated_by: spec/check_vectors.py --emit`. No such flag
exists, and no tool generates the file. The claim is removed.

## Correction 5 — the drift gap comes from dose coupling, not geometry

**2026-09-18.** The G-D note in section 5 attributes the gap between our drift and
the study's to "different geometry". The measurement says otherwise.

Our reduced arm gives a closed drift of −1.4290 at 2000 sites. Re-run at the
study's own 13978 sites it gives −1.4163, a 0.9% change. The population size is
not the cause.

The cause is how tightly dose follows expression. With no transport code, dose is
exactly proportional to expression. Real transport deposits energy over a finite
range, which blurs that link. Writing the coupling as
`w_i ∝ rho * e_i + (1 - rho) * mean(e)`:

| rho | closed drift |
|---|---|
| 1.0, as committed | −1.4290 |
| 0.6 | −0.9932 |
| 0.3 | −0.6144 |
| 0.2 | −0.5181 |
| 0.0 | −0.1362 |

The study's −0.4838 falls at rho between 0.15 and 0.2. G-D is unaffected, because
it tests sign and ordering, not magnitude.

## Correction 6 — Correction 3 claimed more than the uniform-dose form can show

**2026-09-18.** Correction 3 says the uniform-dose G-S form makes "a biased refill
visible". It does not.

With every cell at `e = 1.0` and `sigma_div = 0`, every daughter is exactly 1.0,
whichever parent refills the site. A refill that chooses its parents with any
bias therefore leaves the mean at exactly 1.0, and the gate stays green.

The uniform-dose form still has value. It makes death and refill execute, which
catches a crash, a NaN or a corrupted state. It does not catch a biased refill.

The power to catch one lives in the `--controls` cases, which make `e`
non-uniform. The rest of Correction 3 stands: the zero-dose form alone runs no
death and no refill at all.

## Correction 7 — gate G-B could not fail

**2026-09-18.** G-B in `spec/check_selection.py` was arithmetic alone. It toggled
a local flag inside its own loop and never called the cycle code.

So it proved only that the LQ formula is non-linear, which is always true. No
path it claimed to guard could bank dose. A port reproduced it faithfully and
shipped the same decoration.

**What changed.** G-B now drives the real death sweep in `run_arm` under uniform
dose. With dose consumed at the check, the mean SF applied equals the one-cycle
SF every cycle. Banking pulls cycles two to four below it:

```
dose consumed : 6.964060e-02  6.964060e-02  6.964060e-02  6.964060e-02
dose banked   : 6.964060e-02  6.543330e-02  6.576581e-02  6.429799e-02
```

Removing the reset from the default path fails G-B, which is the check that the
old version could never make.

## Correction 8 — the checker and section 4 disagreed at the survival boundary

**2026-09-18.** Section 4 says a cell divides when `RNG <= SF_i` and dies when
`RNG > SF_i`. `spec/check_selection.py` killed on `RNG >= SF` at three sites.

The two rules differ only when the draw equals SF exactly. At `SF = 0` the spec
lets a draw of exactly 0.0 survive, and the checker killed it.

The checker now follows section 4. Its output is byte-identical before and after,
so no published number moved.

One consequence follows from section 4 itself. At 1000 Gy, SF underflows to
exactly 0.0, so a draw of exactly 0.0 survives with probability 2^-53 per cell.
G-Q's "every exposed cell dies" is therefore exact for its pinned seed, not a
theorem.

## Correction 9 — the 50x share cap refuses nothing below 51 cells

**2026-09-18.** Amendment A4 sets a 50x cap on the single-cell share of the
administered activity, with no minimum population.

The share is `max(e) / mean(e)`. One cell with all the expression gives exactly
`n`, the living population. So the share never exceeds `n`, and at 50 cells or
fewer the cap can never refuse. A port running 14 to 48 cells found G-H green for
every uptake rule.

**The rule, corrected.** The cap applies only when more than 50 cells are alive.
Below that, report the share and do not treat the cap as a gate.

`spec/check_selection.py` now refuses to report G-H as a pass when the population
cannot exceed the cap. A control fixes the boundary: 50 cells unreachable, 51
reachable.

---

# Amendments

A red-team read of the circulated spec on 2026-09-17 raised four ways an exact
implementation could be wrong while passing every gate in section 5. Three
reproduce. One does not, and the reason it does not is recorded here as well.

Each amendment changes a rule, so it is separate from the corrections above,
which only fixed descriptions of the reference.

## A1 — the survival check is forced at the end of the cycle

**Confirmed.** Section 4 draws the survival coin "at the first post-exposure
division attempt". In a crowded lattice, `lambda_V` and `lambda_S` can hold a
cell below its division volume for the whole window. The cell then never
attempts division, the coin is never drawn, and it absorbs any dose while
`cell_state` stays 1.

The committed reference does not have the hole. Its death sweep is
unconditional and runs before the division loop
(`study3_closed_loop/study3_run.py:142-149`).

**Rule.** Every cell exposed to dose is drawn exactly once per cycle. The draw
happens at the cell's first division attempt, or at the cycle's final Monte
Carlo step when no attempt occurred. Gate `G-Q` in `spec/check_selection.py` arrests 30% of the
population below division volume at 1000 Gy and requires all 2000 deaths.

## A2 — `accumulated_dose` is consumed by the check, and daughters start at zero

**Confirmed.** Section 1 reads "Set once per cycle from the
transport/dose-accumulation step. Not decremented." A reader can take the field
name and "not decremented" together as a running total across cycles.

The running total breaks the model. `G(T)` is derived for one continuous exposure of
duration `T`. PRRT cycles are about eight weeks apart and repair completes
between them, so the summed dose may not be squared.

Measured cost of the wrong reading, at the pinned alpha and beta:

| Two cycles | Correct `SF(a)*SF(b)` | Banked `SF(a+b)` | Ratio |
|---|---|---|---|
| 2 + 2 Gy | 3.748788e-01 | 3.670324e-01 | 0.979 |
| 5 + 5 Gy | 7.948366e-02 | 6.964060e-02 | 0.876 |
| 10 + 10 Gy | 4.849813e-03 | 2.858008e-03 | 0.589 |

The reference computes dose fresh each cycle and never banks
(`study3_run.py:226-232`).

**Rule.** `accumulated_dose` is reset to 0.0 by the survival check, whether the
cell lived or died. Daughters are born with `accumulated_dose = 0.0`. The LQ
formula sees only the current cycle's dose. Gate `G-B` requires the per-cycle
`SF` to be identical across four equal-dose cycles.

## A3 — a dying cell below two sites is deleted outright

**Plausible, and not covered by a runnable gate here.** Semantics B sets
`V_target = 0` and lets Metropolis dynamics resorb the cell. Setting the target
to zero does not guarantee removal: a favourable adhesion term can make the
shrinking copy attempts unprofitable, and the cell can persist at one or two
sites as a barrier that still blocks refill.

Study 3 has no Potts Hamiltonian, so nothing in this repository can reproduce
the failure. The upstream CPM is the evidence instead: `biofilms_potts.jl:647`
culls only at `cell.volume <= 0`, with no threshold above zero.

**Rule.** Under Semantics B, a cell in the dying state whose actual volume falls
below 2 sites is deleted from the lattice immediately, and its sites are freed.
Log the deletion; a resorption that never completes is a result, not a detail.

## A4 — hoarding is instrumented, not clamped

**Does not reproduce.** The red-team reading is that multiplicative drift is an
unbounded walk in log space, so a cell climbs to `e = 1e10` and takes nearly all
the administered activity.

Measured over 4 cycles at the pinned `sigma_div = 0.10`, expression walks
**down**, not up:

| Seed | max e after 4 cycles | min e | deepest lineage |
|---|---|---|---|
| 20261001 | 0.3827 | 0.021238 | 17 |
| 20261002 | 0.5302 | 0.029513 | 16 |
| 20261003 | 0.4277 | 0.048845 | 18 |

The initial maximum is 5.7942, so the largest value falls by more than tenfold.
`e = 1e10` needs `ln e = 23`, about 52,900 divisions of 0.10 drift along one
lineage; the deepest lineage stops at 18.

The loop suppresses the walk it is accused of amplifying. A high-`e` cell takes
more activity, receives more dose, and dies. Downward drift is the model's
result, not a defect, and an invented `e_max` would be a declared parameter with
no measured provenance.

The residual risk is real but different: an uptake rule that concentrates. Max
single-cell share, in multiples of the uniform share:

| Uptake rule | Worst share |
|---|---|
| proportional to `e`, as specified | 5.8x |
| proportional to `e**4`, a concentrating bug | 132.1x |

**Rule.** No clamp on `expression_e`. Instead, report the maximum single-cell
share of the administered activity every cycle, and refuse above 50x the
uniform share. Gate `G-H` enforces it, and its control is the `e**4` rule.

*(Correction 9 limits this rule to populations above 50 cells.)*

---

# Answers to the four questions in issue #1

## Answer 1 — `SIGMA_DIV = 0.10` and the volume-gated trigger are both confirmed

**The pin is correct.** `sigma_div = 0.10` at `study3_closed_loop/study3_run.py:47`,
recorded in `study3_closed_loop/verdict.json` under `declared`, and stated in
the repository README among the declared parameters. Section 1 of this spec now
names it.

**The volume trigger is not an adaptation. It is the specification.** The
preprint, which the README names as the citation of record, describes exactly
it. `preprint/closing_the_loop_prRT_spatial_cpm.tex:445-446`:

> `$H_{\text{growth}}$` drives target volume toward doubling over the
> clone-specific doubling time, at which point division is attempted and
> Eq.~\ref{eq:sf} is applied.

and at line 403, dying cells are "cleared by **volume resorption**". Those two
rules are Semantics B. The `P_STEP` Poisson clock in `study3_run.py:51` is the surrogate
study 3 used because it had no Potts dynamics, authorized as a declared
boundary in `study3_closed_loop/PROTOCOL.md:185-187`.

**One coupling to be careful with.** In the reference, "death at first
division" is a label, not a coupling. The death draw is an unconditional
instantaneous sweep over every occupied cell at the start of the window, before
the division loop runs (`study3_run.py:142-149`). It does not wait for a cell
to reach a division.

If you gate the death draw on a cell actually reaching `volume >= 2*V_target`,
then a cell that never grows that far inside the window survives by
default. The extinction gate `K3` requires 13978 deaths of 13978 at 1000
Gy/cycle and would break. Keep the death draw unconditional, and let the volume
gate schedule only the refill.

**One gate to add.** A volume trigger sets the division rate from the Potts
growth dynamics, not from a clock with mean 58.4 h. Report your realized
population doubling time and require it inside a declared band around 58.4 h.
Without that, the volume gate can silently change how many times SF is applied
per window, and therefore the kill per cycle, with no gate reporting it.

## Answer 2 — heterogeneous, from the hashed geometry file; not uniform e = 1

`e_i(T0)` is the T0 activity map of the pinned geometry object, not a uniform
value. `study3_run.py:99`:

```python
E_T0 = gdat["activity"].astype(np.float64).ravel(order="C").copy()
```

- File: `data/geometry_pitch40.npz` in this repository, key `activity`, float64,
  shape (25, 25, 25).
- sha256 `a6883b84a577dd5db888c9d4292cf68897a68d840b64a427061bc70a2d60a170`,
  asserted before use at `study3_run.py:89-92`.
- Generated as lognormal with `mu = 0.0`, `sigma = 0.6`, seed 20260917, then
  rescaled so the arithmetic mean over occupied voxels is exactly 1.0.
- Measured on the committed array: 13978 viable voxels, mean 1.0, min 0.0863,
  max 8.5914, CV 0.658197, mean `ln e` -0.17847, sd `ln e` 0.59654.

Only relative values matter: drift is multiplicative and uptake is
renormalized each cycle.

Uniform `e = 1` appears in exactly one place, the `K2` null arm
(`study3_run.py:182`), which is the G-S control. Use it there and nowhere else.

## Answer 3 — the current `.bin` plus sidecar, unchanged

**Transport app → response layer.** `<out>.bin` is float64, `n^3` values,
C-order, MeV per voxel, no header. Voxel index is `(i*n + j)*n + k` with i along
x. `<out>.json` is the sidecar.

The response layer reads exactly one key from the sidecar,
`decays_simulated`, and uses it once as a divisor (`study3_run.py:111`):

```python
E = np.fromfile(out + ".bin", dtype=np.float64) / meta["decays_simulated"]
```

The divisor converts MeV/voxel-total into MeV per decay per voxel. `n` and `pitch_um`
are present in the sidecar but ignored on read. The array size is checked
against the geometry instead, with `assert E.size == NV`.

Gy conversion is a second, separate step. `kphys` is computed once at the first
transport of the study so that mean viable dose is 10 Gy, then frozen across
every arm and cycle. Recorded value:
`kphys_Gy_per_MeV_per_decay = 1528990.2432354644` in `study3_closed_loop/verdict.json`.

**Response layer → transport app**, for completeness. `<tag>.bin` is int32
`cell_id[n^3]` C-order immediately followed by float64 `activity[n^3]` C-order,
one file, no header; sidecar `<tag>_meta.json` is `{"n": 25, "pitch_um": 40.0}`.
For n=25 that is 187500 bytes. The transport app reads `cell_id` and never uses
it; only `activity` matters, as the sampling CDF.

No new format. Your section-7 assumption is correct: dose arrives
pre-accumulated per cycle, so no MCS-to-hours mapping enters the module.

## Answer 4 — statistical agreement only, and the test should be declared

Confirmed. Nothing in the specification requires bit-matched normal draws, and
nothing downstream can see the difference:

- Division contributes zero expected `ln e` drift. The daughter's value is
  `e_parent * exp(sigma_div * N(0,1))` and the parent is left untouched, so
  `E[ln e_daughter] = E[ln e_parent]`. Every directional drift in G-D comes
  from the death sweep, not from the normals.
- `G-S` runs at `sigma_div = 0`, where `exp(0 * z) == 1.0` bit-exactly for any
  `z`. The RNG cannot reach that gate.
- `G-D` is a sign-and-ordering criterion, not a value criterion, and the
  measured margin is wide: -0.4838 closed against -0.1164 for the worst open
  replicate, contrast +0.2952.

Box-Muller against numpy's Ziggurat is therefore fine. One request: declare the
statistical test rather than asserting that agreement is sufficient. State the
seed count, the statistic and the band. "Statistically sufficient" with no
stated band is not a gate, because nothing can fail it.

Expect no value in `verdict.json` to reproduce bit-for-bit in any case. Changing
the trigger changes RNG stream consumption, so every recorded number moves. The
verdicts are what should survive.
