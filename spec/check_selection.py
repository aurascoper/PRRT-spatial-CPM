#!/usr/bin/env python3
"""Run gates G-S, G-D and the amendment gates on a reduced population.

G-N, G-O, G-P and G-C are closed-form and live in spec/check_vectors.py. The
two gates here need a population that dies and refills over several cycles, so
they run a toy arm instead: same mechanism as study3_closed_loop/study3_run.py,
smaller lattice, no transport.

The absolute drift numbers are not the study's and are not meant to be. The
spec says so for G-D: the sign and the closed > open ordering are the physics.

    python3 spec/check_selection.py             # the two gates, exit 0 or 1
    python3 spec/check_selection.py --controls  # also prove each gate can fail

Stdlib only, Python 3.12 or later. Runs in about a second.
"""
import argparse
import math
import random
import sys

if sys.version_info < (3, 12):  # not an assert: python -O strips those
    sys.exit("Python 3.12+ needed: random.binomialvariate")

# ---- declared constants, spec section 2 and study3_run.py:44-52 -------------
ALPHA, BETA = 0.24, 0.06
T_REP, T_ACT = 1.5, 96.0
MU_REP = math.log(2.0) / T_REP
SIGMA_DIV = 0.10
T_D = 58.4
DT_DAY = 0.25
P_STEP = 1.0 - math.exp(-(math.log(2.0) / T_D) * 24.0 * DT_DAY)
N_STEPS = 224
D_MEAN_TARGET = 10.0
A_ADMIN = 1.0

# ---- reduced-population declarations (not the study's) ---------------------
CAP = 2000          # study 3 uses 13978 viable sites
EXPR_SIGMA = 0.6    # same lognormal spread as data/geometry_pitch40.npz
N_CYCLES = 4
SEEDS = (20261001, 20261002, 20261003, 20261004, 20261005)
CONTRAST_MARGIN = 0.02
ARREST_FRAC = 0.30      # share of cells held below division volume, gate G-Q
HOARD_MULT = 50.0       # max single-cell activity share, in multiples of uniform
LETHAL_GY = 1000.0


def G_of_T(T_h):
    x = MU_REP * T_h
    if x < 1e-3:
        return 1.0 - x / 3.0 + x * x / 12.0
    return 2.0 * (x + math.expm1(-x)) / (x * x)


G96 = G_of_T(T_ACT)


def initial_expression(seed=20260917):
    """Lognormal, then rescaled so the mean over sites is exactly 1.0."""
    r = random.Random(seed)
    e = [r.lognormvariate(0.0, EXPR_SIGMA) for _ in range(CAP)]
    m = sum(e) / CAP
    return [v / m for v in e]


E_T0 = initial_expression()
# kphys is set once so mean viable dose is D_MEAN_TARGET, then frozen, as in
# study3_run.py:226-231.
KPHYS = D_MEAN_TARGET / (A_ADMIN / CAP)


def run_arm(seed, closed, sigma_div=SIGMA_DIV, uniform_e=False,
            dose="hetero", inherit=None, uptake=None,
            bank_dose=False, sf_log=None):
    """One arm of n_cycles. Returns (mean ln e per cycle start, deaths per cycle).

    closed     True re-normalizes activity over the living each cycle.
               False freezes the cycle-1 voxel map and re-applies it blind, so a
               dead site keeps its declared activity (study3_run.py:20-23).
    dose       "hetero", "uniform" (every site D_MEAN_TARGET) or "zero".
    inherit    daughter rule, for the negative controls.
    uptake     weight rule, for the negative controls.
    bank_dose  the A2 defect, for gate G-B's control: carry each site's dose
               into the next cycle instead of consuming it at the check.
    sf_log     if a list, receives the mean SF applied to the living each cycle.
    """
    rng = random.Random(seed)
    inherit = inherit or (lambda ep, s, r: ep * math.exp(s * r.gauss(0.0, 1.0)))
    uptake = uptake or (lambda ev: ev)

    e = [1.0] * CAP if uniform_e else list(E_T0)
    acc = [0.0] * CAP      # dose the survival check reads; A2 says it is consumed
    occupied = list(range(CAP))
    free = []
    frozen_w = None
    ln_means, death_counts = [], []

    for _ in range(N_CYCLES):
        if not occupied:
            break
        ln_means.append(sum(math.log(e[i]) for i in occupied) / len(occupied))

        if dose == "zero":
            D = [0.0] * CAP
        elif dose == "uniform":
            D = [D_MEAN_TARGET] * CAP
        else:
            w = [0.0] * CAP
            for i in occupied:
                w[i] = uptake(e[i])
            tot = sum(w)
            w = [v / tot * A_ADMIN for v in w] if tot > 0 else w
            if frozen_w is None:
                frozen_w = list(w)
            D = [KPHYS * v for v in (w if closed else frozen_w)]

        # A2: the check reads this cycle's dose only. Banking it instead is the
        # defect gate G-B exists to refuse, so it is reachable here on purpose.
        for i in occupied:
            acc[i] = acc[i] + D[i] if bank_dose else D[i]

        # death sweep: unconditional, before any division (study3_run.py:142-149)
        died = []
        sf_sum = 0.0
        for i in list(occupied):
            sf = math.exp(-(ALPHA * acc[i] + BETA * G96 * acc[i] ** 2))
            sf_sum += sf
            if rng.random() > sf:          # spec §4: dies iff RNG > SF
                died.append(i)
        if sf_log is not None:
            sf_log.append(sf_sum / len(occupied))
        dead = set(died)
        occupied = [i for i in occupied if i not in dead]
        free.extend(died)
        death_counts.append(len(died))

        # refill to capacity
        for _ in range(N_STEPS):
            if not free or not occupied:
                break
            k = min(rng.binomialvariate(len(occupied), P_STEP), len(free))
            if k == 0:
                continue
            parents = rng.sample(occupied, k)
            for p in parents:
                j = rng.randrange(len(free))
                free[j], free[-1] = free[-1], free[j]
                t = free.pop()
                e[t] = inherit(e[p], sigma_div, rng)
                acc[t] = 0.0       # A2: a daughter is born with no dose
                occupied.append(t)
    return ln_means, death_counts


# ---- gate G-S ---------------------------------------------------------------

def gate_S(inherit=None):
    """Uniform e, zero drift: mean e must stay exactly 1.0, with and without dose.

    The spec says zero dose. The committed K2 arm uses uniform dose instead
    (study3_run.py:253-256), so death and refill actually execute. Both are
    required here.

    WHAT THE UNIFORM-DOSE FORM CANNOT SEE. With every cell at e = 1.0 and
    sigma_div = 0, every daughter is exactly 1.0 whichever parent refills the
    site. So a refill that picks its parents with any bias still leaves the mean
    at exactly 1.0, and this main assertion stays green. The uniform-dose form
    catches a crash, a NaN or a corrupted state during death and refill. It does
    not catch a biased refill. That discriminating power lives only in the
    --controls cases, which make e non-uniform. An earlier docstring, and the
    spec's Correction 3, claimed the opposite.
    """
    detail = []
    for label, dose in (("zero dose", "zero"), ("uniform dose", "uniform")):
        ln_means, deaths = run_arm(20261201, closed=True, sigma_div=0.0,
                                   uniform_e=True, dose=dose, inherit=inherit)
        ok = all(m == 0.0 for m in ln_means)     # ln(1.0) is exactly 0.0
        detail.append(f"{label}: deaths/cycle {deaths}, mean ln e {ln_means}")
        if not ok:
            return False, detail
    return True, detail


# ---- gate G-D ---------------------------------------------------------------

def gate_D(uptake=None, open_is_closed=False):
    """Heterogeneous dose: mean expression falls, and closed falls further.

    The two arms share a dynamics seed, so the comparison is paired. Cycle 1 is
    identical in both arms anyway, because the open arm freezes the cycle-1 map.
    Pairing removes the seed-to-seed spread that an unpaired mean cannot see:
    an earlier version of this gate compared unpaired means against a fixed
    0.02 margin and accepted its own negative control.
    """
    closed_d, open_d, per_seed = [], [], []
    for s in SEEDS:
        c, _ = run_arm(s, closed=True, uptake=uptake)
        o, _ = run_arm(s, closed=open_is_closed, uptake=uptake)
        dc, do = c[-1] - c[0], o[-1] - o[0]
        closed_d.append(dc)
        open_d.append(do)
        per_seed.append(abs(dc) - abs(do))
    mc = sum(closed_d) / len(closed_d)
    mo = sum(open_d) / len(open_d)
    contrast = sum(per_seed) / len(per_seed)
    ok = (mc < 0.0 and mo < 0.0 and contrast > CONTRAST_MARGIN
          and all(v > 0.0 for v in per_seed))
    return ok, [f"closed mean drift {mc:+.4f}", f"open mean drift {mo:+.4f}",
                f"paired contrast {contrast:+.4f} (needs > {CONTRAST_MARGIN})",
                f"per-seed contrast {[round(v, 4) for v in per_seed]} "
                f"(every one must be > 0)"]


# ---- gate G-Q: arrested cells are not immortal ------------------------------

def gate_Q(draw_only_on_division=False):
    """A lethal dose must kill every exposed cell, including arrested ones.

    A cell held below its division volume by crowding never attempts division.
    If the survival draw is gated on that attempt, the cell absorbs any dose and
    stays viable. The committed reference does not have the hole: its death
    sweep is unconditional (study3_run.py:142-149).
    """
    rng = random.Random(20261301)
    n_arrested = int(CAP * ARREST_FRAC)
    arrested = set(range(n_arrested))
    sf = math.exp(-(ALPHA * LETHAL_GY + BETA * G96 * LETHAL_GY ** 2))
    deaths = 0
    for i in range(CAP):
        if draw_only_on_division and i in arrested:
            continue                      # never attempts division, never drawn
        if rng.random() > sf:              # spec §4: dies iff RNG > SF
            deaths += 1
    ok = deaths == CAP
    return ok, [f"{LETHAL_GY:g} Gy, SF = {sf!r}",
                f"{n_arrested} of {CAP} cells arrested below division volume",
                f"deaths {deaths} of {CAP} (every exposed cell must die)"]


# ---- gate G-B: dose is consumed by the check, never banked ------------------

def gate_B(bank_dose=False):
    """SF at cycle k must use cycle k's dose, not the sum over earlier cycles.

    G(T) is derived for one continuous exposure of duration T. PRRT cycles are
    weeks apart and repair completes between them, so squaring a summed dose is
    not the same model.

    THIS GATE USED TO BE ARITHMETIC ALONE, and could not fail. It looped over a
    local list, toggling a local `accumulate` flag, and never called run_arm. So
    it showed only that the LQ formula is non-linear, which is always true, and
    no code path it guarded could bank dose. A port that reproduced it faithfully
    shipped the same decoration, which is how the defect was found.

    It now drives run_arm's real death sweep under uniform dose. With dose
    consumed at the check, every living cell sees exactly D_MEAN_TARGET each
    cycle, so the mean SF applied is the same number every cycle. Banking makes
    each survivor carry its old dose forward, and the mean SF falls.
    """
    sf_log = []
    run_arm(20261401, closed=True, dose="uniform", bank_dose=bank_dose,
            sf_log=sf_log)
    want = math.exp(-(ALPHA * D_MEAN_TARGET + BETA * G96 * D_MEAN_TARGET ** 2))
    ok = len(sf_log) == N_CYCLES and all(abs(s - want) < 1e-12 for s in sf_log)
    return ok, [f"mean SF applied per cycle {[f'{s:.6e}' for s in sf_log]}",
                f"required every cycle: {want:.6e}, the SF of one cycle's dose"]


# ---- gate G-H: no cell takes an implausible share of the activity ----------

def share_cap_is_reachable(n, cap=HOARD_MULT):
    """Whether a share cap can ever refuse in a population of n living cells.

    The share is max(e) / mean(e), and one cell holding all the expression gives
    exactly n. So the share can never exceed n, and a cap at or above n refuses
    nothing: gate G-H passes vacuously for any uptake rule at all. A port running
    14 to 48 cells hit exactly that, and the spec never said so.
    """
    return n > cap


def gate_H(power=1.0):
    """Uptake is renormalized, so one cell hoarding it starves every other.

    The threshold is measured, not declared from biology: uptake proportional
    to e peaks at 5.8x uniform, while uptake proportional to e**4 reaches 132x.
    """
    if not share_cap_is_reachable(CAP):
        return False, [f"VACUOUS: the share can never exceed the population "
                       f"{CAP}, so a {HOARD_MULT:g}x cap refuses nothing"]
    uniform = 1.0 / CAP
    worst, per_seed = 0.0, []
    for s in SEEDS:
        rng = random.Random(s)
        e = list(E_T0)
        occupied, free = list(range(CAP)), []
        top = 0.0
        for _ in range(N_CYCLES):
            if not occupied:
                break
            w = [0.0] * CAP
            for i in occupied:
                w[i] = e[i] ** power
            tot = sum(w)
            w = [v / tot for v in w]
            top = max(top, max(w))
            D = [KPHYS * v for v in w]
            died = [i for i in occupied
                    if rng.random() > math.exp(   # spec §4: dies iff RNG > SF
                        -(ALPHA * D[i] + BETA * G96 * D[i] ** 2))]
            dead = set(died)
            occupied = [i for i in occupied if i not in dead]
            free.extend(died)
            for _ in range(N_STEPS):
                if not free or not occupied:
                    break
                k = min(rng.binomialvariate(len(occupied), P_STEP), len(free))
                if k == 0:
                    continue
                for pa in rng.sample(occupied, k):
                    j = rng.randrange(len(free))
                    free[j], free[-1] = free[-1], free[j]
                    tg = free.pop()
                    e[tg] = e[pa] * math.exp(SIGMA_DIV * rng.gauss(0.0, 1.0))
                    occupied.append(tg)
        per_seed.append(top / uniform)
        worst = max(worst, top / uniform)
    ok = worst < HOARD_MULT
    return ok, [f"max single-cell share, in multiples of uniform: "
                f"{[round(v, 1) for v in per_seed]}",
                f"worst {worst:.1f}x (limit {HOARD_MULT:g}x)"]


# ---- controls ---------------------------------------------------------------

def _drift_ignores_sigma(ep, s, r):
    """Hardcodes the drift, so sigma_div = 0 does not silence it. G-S must bite."""
    return ep * math.exp(SIGMA_DIV * r.gauss(0.0, 1.0))


def _inherit_offset(ep, s, r):
    """Additive offset on the daughter. G-S must bite."""
    return ep * math.exp(s * r.gauss(0.0, 1.0)) + 1e-3


def _uptake_inverse(ev):
    """Activity inversely proportional to expression. G-D must bite on sign."""
    return 1.0 / ev


def run_controls():
    cases = [
        ("G-Q", lambda: gate_Q(draw_only_on_division=True),
         "survival drawn only when a cell attempts division"),
        ("G-B", lambda: gate_B(bank_dose=True),
         "dose banked across cycles on the real sweep"),
        ("G-H", lambda: gate_H(power=4.0),
         "uptake proportional to e**4"),
        ("G-S", lambda: gate_S(inherit=_drift_ignores_sigma),
         "drift hardcoded, ignores sigma_div = 0"),
        ("G-S", lambda: gate_S(inherit=_inherit_offset),
         "daughter gets an additive offset"),
        ("G-D", lambda: gate_D(uptake=_uptake_inverse),
         "uptake inversely proportional to e"),
        ("G-D", lambda: gate_D(open_is_closed=True),
         "open arm re-normalized like the closed arm"),
    ]
    fails = []
    print("negative controls (each gate must REFUSE its broken input):")
    for name, fn, why in cases:
        ok, _ = fn()
        print(f"  {'refused ' if not ok else 'ACCEPTED'} {name:<4} {why}")
        if ok:
            fails.append(f"{name} accepted: {why}")
    return fails + check_share_cap_control()


def check_share_cap_control():
    """The vacuity guard on G-H is itself a guard, so it needs a control.

    The boundary is the whole content: at n = cap the most extreme possible
    share equals the cap and cannot exceed it, so n = 50 must read unreachable
    and n = 51 reachable.
    """
    cases = [(14, False), (40, False), (50, False), (51, True), (CAP, True)]
    fails = []
    print("share-cap reachability (a cap at or above the population refuses nothing):")
    for n, want in cases:
        got = share_cap_is_reachable(n)
        ok = got == want
        print(f"  {'ok  ' if ok else 'FAIL'} n={n:<5} reachable={got}  (want {want})")
        if not ok:
            fails.append(f"share-cap reachability wrong at n={n}")
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", action="store_true",
                    help="also prove every gate can fail")
    args = ap.parse_args()

    print(f"reduced population: {CAP} sites, {N_CYCLES} cycles, "
          f"{len(SEEDS)} seeds, G(96 h) = {G96:.6f}")
    fails = []
    for name, fn in (("G-S", gate_S), ("G-D", gate_D),
                     ("G-Q", gate_Q), ("G-B", gate_B), ("G-H", gate_H)):
        ok, detail = fn()
        print(f"  {'PASS' if ok else 'FAIL'} {name}")
        for d in detail:
            print(f"       {d}")
        if not ok:
            fails.append(name)
    if args.controls:
        fails += run_controls()

    print()
    if fails:
        print(f"REFUSED: {len(fails)} check(s) failed")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
