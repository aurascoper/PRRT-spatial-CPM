# What the `response.growth_survival` PR must contain

Acceptance criteria for the pull request that implements
`spec/response_growth_survival.md` in a Cellular Potts port. Written against
the Odin port at `ffinkdevs/Biofilms`, tag `odin-ensemble-256`, commit
`4788dcd1`, so the paths and hazards below are that port's actual ones. Reviewed
at `f3bf6d6` on branch `odin-port`; the tag still points at `4788dcd1`.

The spec says what to build. The criteria below say what the pull request has
to show before anyone can tell whether it was built correctly.

## 1. Declarations the PR description must state

A reviewer cannot infer any of these from a diff. State each one in the PR
body, not in a comment thread.

| Declaration | Allowed values |
|---|---|
| Death and refill semantics | `A` (clear immediately, well-mixed refill) or `B` (`V_target = 0`, Metropolis resorption) |
| Division trigger | volume-gated, or time-stochastic with the declared `T_D` |
| RNG stream for survival draws | which generator, and in what order cells are drawn |
| RNG stream for drift normals | which generator, and which normal algorithm |
| `sigma_div` | the pinned `0.10`, or a declared deviation with its reason |
| Amendments implemented | A1, A2, A3, A4 from the spec, each yes or no |

Mixing semantics A and B silently fails review. Declaring B and implementing
A's refill also fails review, and the gate ladder cannot catch it, so the
declaration is the only record.

## 2. Four hazards specific to this port

Each was read from the port at `4788dcd1`. None is hypothetical. Line numbers in
the hazard text are at that commit. Each hazard ends with a status line from the
review of `f3bf6d6`. The port's 26 tests were run there, then re-run against a
planted defect wherever a status line says so.

### H1 — three memory layouts, three places to add every field

`odin/cpm/sim.odin:10` gives `Sim` one of `aos`, `soa` or `aosoa`, and every
accessor branches on `s.layout`. `cell_alive`, `cell_kill` and `cell_set_new`
each switch three ways.

The spec adds three per-cell fields: `accumulated_dose`, `cell_state` and
`expression_e`. Add each to all three layouts. A field added to `AoS` alone
makes the port's answer depend on a build flag.

`odin/tests/layouts_test.odin` already asserts that the three layouts agree on
volume and melanin after 5 MCS. Extend `test_layout_equivalence` to run one
full dose cycle and compare the new fields, or the existing test will pass over
a port that disagrees with itself.

**Status at `f3bf6d6`: partly closed.** `odin/cpm/layouts.odin` adds `dose`,
`state` and `expr` to all three layouts, with a `grow_registry_*` for each.
`test_layout_equivalence` runs one 5 Gy cycle and compares the three fields per
cell. Two of the three comparisons cannot fail where they stand. Expression is
never initialised or inherited in that test, so it compares 1.0 with 1.0. Dose
is compared after the A2 reset, so it compares 0.0 with 0.0.

Run against the port: an `SoA` `cell_set_expr` that writes nothing passes all 26
tests, and so does an `SoA` `cell_expr` that reads slot 0. A dose accessor that
returns 0 is caught, through `state`. The cycle divides no cell, so
`grow_registry_soa` and `grow_registry_aosoa` never run in any test.

Three lines before the cycle close it: `init_expression_lognormal` on each sim,
`assign_dose_scaled_e` in place of the uniform dose, and `v_target = 2` so cells
divide. With those, both expression defects fail the test. The port as committed
passes it, so its three layouts agree through death, division and registry
growth.

### H2 — the survival draws must not consume the Julia stream

`Sim` declares two generators (`odin/cpm/sim.odin:15-17`):

```odin
rng:      Rng,        // splitmix64 fast stream
jr:       Julia_Rng,  // Julia 1.12 MersenneTwister stream (opt-in)
```

The 256-seed ensemble at `odin/compare/ensemble_n40_p6_42-297_odin.csv` matches
the Julia reference pathwise because `jr` replays Julia's draw sequence exactly.
Any new draw from `jr` shifts that sequence and every parity result dies with
it.

**Rule.** Survival uniforms and drift normals come from `s.rng`. If the PR needs
them from `jr` instead, it must re-pin the parity fixtures in the same PR and
say so in the title.

**Status at `f3bf6d6`: closed, by declaration and by reading.** The commit message
declares both streams as `s.rng`, never `jr`. `odin/cpm/response.odin` draws at
lines 332, 360 and 361, all three from `s.rng`, and names `jr` in comments only.
No test was run against a defect here.

### H3 — there is no divide primitive yet

`cell_kill` exists at `odin/cpm/sim.odin:229`. There is no `cell_divide`. The
PR builds it.

A new cell must update `next_id`, `n_slots`, and the per-layout `n_alive`
counters. `SoA` and `AoSoA` maintain `n_alive` by hand, and `AoS` does not
maintain it at all, so a divide that increments one and not the others makes
`cell_alive` and the population count disagree.

**Status at `f3bf6d6`: closed in code, exercised under `AoS` only.** `cell_divide`
is at `odin/cpm/sim.odin:390-462`. Daughters take their counters through
`cell_set_new`, and ids are never reused. Every response test builds an `AoS`
sim, and the layout test divides nothing, so no committed test divides a cell
under `SoA` or `AoSoA`. The three lines under H1 cover that as well.

`cell_divide` checks that the parent is alive and checks no state. Amendment A5
in the spec is about the one line that stands in for that check.

### H4 — `cell_kill` does not clear the lattice, and a dead sigma spreads

H4 is the reason A3 matters more in this port than in the reference.

`cell_kill` flips `alive` and touches nothing else. The cell's sites keep its
sigma in `s.arena.lattice`. `species_of` then reports `MEDIUM` for those sites
(`odin/cpm/sim.odin:182-186`), so adhesion treats them as medium.

Today that is harmless, by accident. The reap loop kills only at
`cell_volume(s, id) <= 0`, and volume counts sites, so a cell is killed exactly
when it has no sites left. No dead-sigma site can exist.

**Semantics B breaks the accident.** Amendment A3 culls a dying cell below 2
sites, which kills it while it still occupies one. The site then keeps a dead
sigma, and the copy path propagates it:

```odin
if accept {
    if sig_t > 0 && cell_alive(s, int(sig_t)) { cell_add_volume(s, int(sig_t), -1) }
    if sig_s > 0 && cell_alive(s, int(sig_s)) { cell_add_volume(s, int(sig_s), +1) }
    s.arena.lattice[ti] = sig_s
}
```

A copy whose source is the dead cell's site writes that dead sigma into the
target. The volume update is skipped, because the guard requires `cell_alive`.
The reap loop skips it too, because the guard requires `cell_alive`. So the dead
sigma can occupy any number of sites with nothing tracking it.

**Rule.** The A3 cull writes `0` into `arena.lattice` for every site of the
dying cell before calling `cell_kill`. Assert afterwards that no lattice entry
names a dead cell.

**Status at `f3bf6d6`: closed.** `cull_dying` at `odin/cpm/sim.odin:468-486` zeroes
the sites before `cell_kill`, and `assert_no_dead_sigmas` is at `sim.odin:502`.
With the zeroing removed, `test_resp_A3_cull` fails at
`odin/tests/response_test.odin:562`. `test_resp_H4_no_dead_sigmas` plants the
old kill-without-clear path and requires the invariant to refuse it.

## 3. Files the PR is expected to touch

Matching the port's existing layout, not prescribing a new one.

| Path | Change |
|---|---|
| `odin/cpm/sim.odin` | three per-cell fields in all three layouts, `cell_divide`, the A3 cull |
| `odin/cpm/response.odin` | new: `G(T)`, `SF`, the survival draw, the drift |
| `odin/cpm/params.odin` | `alpha`, `beta`, `T_rep`, `T_active`, `sigma_div`, declared not tuned |
| `odin/tests/response_test.odin` | new: the gate ladder, in the `@(test)` idiom of `julia_exp_test.odin` |
| `odin/tests/layouts_test.odin` | extend `test_layout_equivalence` over the new fields |
| `odin/README.md` | the declarations from section 1 |

## 4. Tests the PR must ship

One test per gate, each with the input that makes it fail. A gate with no such
input reads the same whether it is checking something or nothing.

| Gate | Assertion | The input it must reject |
|---|---|---|
| G-N | `SF(0 Gy) == 1.0` exactly | an `SF` scaled by a constant |
| G-O | `SF(1000 Gy) == 0.0`, never NaN | an `SF` floored at a tiny positive value |
| G-P | `SF(D, 24 h) < SF(D, 96 h)`, and `G(1e-9 h) -> 1` | `G` with no Taylor branch |
| G-C | `12 < 1/G(96 h) < 23` | `T_rep = 15 h` |
| G-S | uniform `e`, `sigma_div = 0`, uniform dose: mean `e` stays `1.0` | drift that ignores `sigma_div = 0` |
| G-D | closed drift below open drift, both negative | uptake inversely proportional to `e` |
| G-Q | every exposed cell dies at 1000 Gy, including arrested ones | the draw gated on a division attempt |
| G-B | per-cycle `SF` identical across four equal-dose cycles | `accumulated_dose` banked across cycles |
| G-H | max single-cell activity share below 50x uniform, with more than 50 cells alive | uptake proportional to `e**4` |
| H1 | the three layouts agree on the new fields after a dose cycle | a field added to `AoS` only |
| H4 | no lattice entry names a dead cell | the cull that calls `cell_kill` without clearing sites |

The vectors for G-N through G-C are in `spec/reference_vectors.json`. The
working forms of every gate are in `spec/check_vectors.py` and
`spec/check_selection.py`, both stdlib-only Python (3.12 or later for the
second), so they can be read as executable pseudocode rather than
reimplemented from prose.

## 5. What makes the PR mergeable

1. Every declaration in section 1 is stated in the PR body.
2. Every gate in section 4 passes, and every rejection input is committed as a
   test that fails when the fix is reverted.
3. The four hazards in section 2 are each addressed or explicitly declined with
   a reason.
4. The existing parity results still hold: the 256-seed ensemble still matches,
   or the PR re-pins the fixtures and says so in the title.
5. `expression_e` is reported per cycle as a distribution, not a mean alone.
   Minimum, maximum, mean of `ln e`, and the maximum single-cell activity share.

## 6. Out of scope for this PR

No bystander signalling, no immune response, no continuous dose-rate
integration inside a cycle, no clonogen expansion past the refill mechanics.

No transport. Dose arrives pre-accumulated per cycle, so no mapping between
Monte Carlo steps and hours enters the module.

Anything beyond the spec needs its own declared provenance in the PR
description.
