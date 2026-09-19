#!/usr/bin/env python3
"""Check an implementation of `response.growth_survival` against the spec.

Asserts every test vector in spec/reference_vectors.json, then runs the four
numeric gates G-N, G-O, G-P and G-C from spec section 5.

Every gate also runs against a deliberately broken implementation that it must
reject. A gate with no such control is decoration: it reads the same whether it
is checking something or nothing. `--controls` prints that ladder.

    python3 spec/check_vectors.py              # vectors + gates, exit 0 or 1
    python3 spec/check_vectors.py --controls   # also prove each gate can fail

Stdlib only. No numpy, no install.
"""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VECTORS = json.load(open(os.path.join(HERE, "reference_vectors.json")))

ALPHA = VECTORS["constants"]["alpha_per_Gy"]
BETA = VECTORS["constants"]["beta_per_Gy2"]
T_REP = VECTORS["constants"]["T_rep_h"]
MU = math.log(2.0) / T_REP


# ---- the reference implementation (spec section 2) --------------------------

def G_of_T(T_h, mu=MU):
    """Lea-Catcheside protraction factor over an exposure of T_h hours."""
    x = mu * T_h
    if x < 1e-3:
        return 1.0 - x / 3.0 + x * x / 12.0        # limit is 1, not 0
    # expm1, never exp(-x). Association as study2_run.py:70, so the bits match it.
    return 2.0 / (x * x) * (x + math.expm1(-x))


def SF(D, G):
    """Surviving fraction after dose D Gy at protraction factor G."""
    return math.exp(-(ALPHA * D + BETA * G * D * D))


# ---- the broken implementations each gate must reject -----------------------

def G_naive(T_h, mu=MU):
    """No Taylor branch. Cancels to 0.0 at short exposure. G-P must catch it."""
    x = mu * T_h
    return 2.0 * (x - 1.0 + math.exp(-x)) / (x * x)


def G_rejected_taylor(T_h, mu=MU):
    """The 0.5*x*(1-x/3) form spec section 2 rejects. G-P must catch it."""
    x = mu * T_h
    return 0.5 * x * (1.0 - x / 3.0) if x < 1e-3 else G_of_T(T_h, mu)


def G_wrong_trep(T_h):
    """T_rep = 15 h instead of 1.5 h. G-C must catch it."""
    return G_of_T(T_h, mu=math.log(2.0) / 15.0)


def SF_floored(D, G):
    """Clamps away the underflow. G-O must catch it."""
    return max(SF(D, G), 1e-300)


def SF_scaled(D, G):
    """Loses the D=0 identity by a hair. G-N must catch it."""
    return SF(D, G) * 0.999


def SF_sign_flip(D, G):
    """+beta instead of -beta. Reverses protraction. G-P must catch it."""
    return math.exp(-ALPHA * D + BETA * G * D * D)


# ---- gates (spec section 5) -------------------------------------------------

def gate_N(sf, g):
    """D = 0 gives SF exactly 1.0."""
    return sf(0.0, g(96.0)) == 1.0


def gate_O(sf, g):
    """D = 1000 Gy underflows to exactly 0.0, never NaN, never negative."""
    v = sf(1000.0, g(96.0))
    return v == 0.0 and not math.isnan(v)


def gate_P(sf, g):
    """Shorter exposure is strictly more lethal, and G tends to 1, not 0."""
    for D in (2.0, 5.0, 10.0, 20.0):
        if not sf(D, g(24.0)) < sf(D, g(96.0)):
            return False
    return abs(g(1e-9) - 1.0) < 1e-6


def gate_C(sf, g):
    """1/G(96 h) falls inside the band the Tamborino PRRT fits imply."""
    band = VECTORS["G_C_band"]
    return band["low"] < 1.0 / g(96.0) < band["high"]


GATES = [("G-N", gate_N), ("G-O", gate_O), ("G-P", gate_P), ("G-C", gate_C)]


# ---- vector assertions ------------------------------------------------------

def close(got, want, sig=6):
    """True when `got` matches `want` to `sig` significant digits.

    Half an ulp at the sig-th significant digit, so `5 * 10**-sig`. An earlier
    version computed `10**(-sig + 1)`, which is a full decade looser than its own
    name: at sig=6 it enforced five digits while every caller and the vectors
    file's `_tolerance` string said six. One published vector, SF(10 Gy, 96 h)
    stored as 0.069641, sat at 5.75e-6 relative and passed only on that slack.
    """
    if want == 0.0 or got == 0.0:
        return got == want
    return abs(got - want) <= abs(want) * 5.0 * 10.0 ** (-sig)


def check_vectors():
    fails = []

    def record(name, got, want, ok):
        mark = "ok  " if ok else "FAIL"
        print(f"  {mark} {name:<34} got {got!r:<24} want {want!r}")
        if not ok:
            fails.append(name)

    print("G(T), spec section 3:")
    for T_str, entry in VECTORS["G_of_T"].items():
        if T_str.startswith("_"):
            continue
        want, sig = entry["value"], entry["sig"]
        got = G_of_T(float(T_str))
        record(f"G({T_str} h) to {sig} sig", round(got, sig + 6), want,
               close(got, want, sig))

    print("branch continuity at x = 1e-3:")
    bc = VECTORS["branch_continuity"]
    x = bc["x"]
    diff = abs((1.0 - x / 3.0 + x * x / 12.0)
               - 2.0 / (x * x) * (x + math.expm1(-x)))
    record("|taylor - standard|", f"{diff:.3e}",
           f"< {bc['tolerance']:.0e}", diff < bc["tolerance"])

    print("SF at G(96 h), spec section 3:")
    g96 = G_of_T(96.0)
    for D_str, want in VECTORS["SF_at_G96"].items():
        got = SF(float(D_str), g96)
        record(f"SF({D_str} Gy)", f"{got:.6e}", f"{want:.6e}", close(got, want))

    print("protraction at 10 Gy:")
    p = VECTORS["protraction_at_10Gy"]
    for T, key in ((24.0, "T_24h"), (96.0, "T_96h")):
        got = SF(10.0, G_of_T(T))
        # FULL PRECISION ON PURPOSE. This line printed round(got, 6), so a value
        # failing its tolerance displayed as identical to its expectation --
        # `got 0.069641  want 0.069641` beside a FAIL. A diagnostic that hides
        # the disagreement it reports is worse than none.
        record(f"SF(10 Gy, T={T:g} h)", got, p[key], close(got, p[key]))

    print("G-C quantity:")
    band = VECTORS["G_C_band"]
    got = 1.0 / g96
    record("1 / G(96 h)", got, band["value"], close(got, band["value"]))
    return fails


def check_gates():
    fails = []
    print("gates on the reference implementation (all must PASS):")
    for name, fn in GATES:
        ok = fn(SF, G_of_T)
        print(f"  {'PASS' if ok else 'FAIL'} {name}")
        if not ok:
            fails.append(name)
    return fails


def check_controls():
    """Each gate against an implementation it must reject."""
    cases = [
        ("G-N", gate_N, SF_scaled, G_of_T, "SF scaled by 0.999"),
        ("G-O", gate_O, SF_floored, G_of_T, "SF floored at 1e-300"),
        ("G-P", gate_P, SF, G_naive, "G with no Taylor branch"),
        ("G-P", gate_P, SF, G_rejected_taylor, "G as 0.5*x*(1-x/3)"),
        ("G-P", gate_P, SF_sign_flip, G_of_T, "beta sign flipped"),
        ("G-C", gate_C, SF, G_wrong_trep, "T_rep = 15 h"),
    ]
    fails = []
    print("negative controls (each gate must REFUSE its broken input):")
    for name, gate, sf, g, why in cases:
        refused = not gate(sf, g)
        print(f"  {'refused' if refused else 'ACCEPTED'} {name:<5} {why}")
        if not refused:
            fails.append(f"{name} accepted: {why}")
    return fails + check_tolerance_control()


def check_tolerance_control():
    """The tolerance is itself a guard, so it gets a control too.

    `close` once enforced five significant digits while claiming six, and the
    only thing that exposed it was an external port. So prove both directions
    here: a value stored at five digits must be refused, the same quantity at
    seven must be accepted, and the old formula must be shown to accept the
    five-digit value -- otherwise nothing proves the tightening changed anything.
    """
    true_sf96 = SF(10.0, G_of_T(96.0))       # 0.0696405992270619
    five_sig = 0.069641                       # what the file once stored
    seven_sig = 6.964060e-02                  # what it stores now

    def old_close(got, want, sig=6):          # the defect, kept as the control
        return abs(got - want) <= abs(want) * 10.0 ** (-sig + 1)

    cases = [
        ("five-digit value refused", not close(true_sf96, five_sig)),
        ("seven-digit value accepted", close(true_sf96, seven_sig)),
        ("old formula accepted the five-digit value", old_close(true_sf96, five_sig)),
    ]
    fails = []
    print("tolerance control (the tolerance must bite, and must differ from the defect):")
    for why, ok in cases:
        print(f"  {'ok  ' if ok else 'FAIL'} {why}")
        if not ok:
            fails.append(f"tolerance: {why}")
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", action="store_true",
                    help="also prove every gate can fail")
    args = ap.parse_args()

    fails = check_vectors() + check_gates()
    if args.controls:
        fails += check_controls()

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
