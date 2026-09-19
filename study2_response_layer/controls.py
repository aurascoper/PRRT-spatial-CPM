#!/usr/bin/env python3
"""Run study2_run.py against three planted defects. Each must end REFUSED-GATED
at the gate named here, or the gate is decoration (issues #11, #13).

    python3 study2_response_layer/controls.py     # exit 0 or 1
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
VERDICT = os.path.join(HERE, "verdict.json")
sha = lambda: hashlib.sha256(open(VERDICT, "rb").read()).hexdigest()
before = sha()
CASES = [("flat", "G6b noise-only null", "flat field with 0.71% MC noise, no heterogeneity"),
         ("pin",  "G7b dose-field hash", "hash of a different real file against the pin"),
         ("unif", "G8 energy bookkeeping", "uniform arm built from the median, not the mean")]
fails = 0
for name, gate, why in CASES:
    r = subprocess.run([sys.executable, os.path.join(HERE, "study2_run.py"), f"--control={name}"],
                       capture_output=True, text=True)
    ok = ("=== VERDICT: REFUSED-GATED ===" in r.stdout and f"GATE {gate}: FAIL" in r.stdout
          and "REFUSED-GATED: verdict.json not written" in r.stdout)
    fails += not ok
    print(f"  {'refused ' if ok else 'ACCEPTED'} {name:5s} {why} -> {gate}")
    if not ok:
        print("      tail:", r.stdout.strip().splitlines()[-2:], r.stderr.strip().splitlines()[-1:])
if sha() != before:
    fails += 1
    print("  ACCEPTED a gated run rewrote verdict.json (the G7b pin would move)")
else:
    print("  refused  gated runs left verdict.json byte-identical")
print("all controls refused" if not fails else f"REFUSED: {fails} control(s) accepted")
sys.exit(1 if fails else 0)
