#!/usr/bin/env python3
"""Run study2_run.py against planted defects. Each must end REFUSED-GATED at the
gate named here, or the gate is decoration (issues #11, #13). G7b is attacked on
the real, unflagged path: a substituted field twice, then with no pin file. That proves a gated run
writes no verdict.json and an absent pin refuses instead of self-pinning.

    python3 study2_response_layer/controls.py     # exit 0 or 1
"""
import hashlib, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
VERDICT = os.path.join(HERE, "verdict.json")
PINS = os.path.join(HERE, "pins.json")
for f in (VERDICT, PINS):
    if not os.path.exists(f):
        sys.exit(f"REFUSED: {f} is missing; restore it from git before running the controls")
sha = lambda f=VERDICT: hashlib.sha256(open(f, "rb").read()).hexdigest()
before, pins_before = sha(), sha(PINS)
CASES = [("flat", "G6b noise-only null", "flat field with 0.71% MC noise, no heterogeneity"),
         ("unif", "G8 energy bookkeeping", "uniform arm built from the median, not the mean")]
fails = 0
for name, gate, why in CASES:
    r = subprocess.run([sys.executable, os.path.join(HERE, "study2_run.py"), f"--control={name}"],
                       capture_output=True, text=True)
    ok = "=== VERDICT: REFUSED-GATED ===" in r.stdout and f"GATE {gate}: FAIL" in r.stdout
    fails += not ok
    print(f"  {'refused ' if ok else 'ACCEPTED'} {name:5s} {why} -> {gate}")
    if not ok:
        print("      tail:", r.stdout.strip().splitlines()[-2:], r.stderr.strip().splitlines()[-1:])
# G7b, first attack: the real gated path, with no --control flag. A garbage field is
# substituted through STUDY2_REF_BIN and the runner runs twice. Run 1 must fail
# G7b; run 2 must fail it again, which only holds if run 1 wrote no verdict.json.
import numpy as np, tempfile
real = os.path.join(HERE, "..", "study1", "runs", "ref_p40_h3.bin")
E = np.fromfile(real, dtype=np.float64)
tmp = tempfile.NamedTemporaryFile(suffix=".bin", delete=False); tmp.close()
(np.random.default_rng(2).uniform(0.1, 10.0, E.size) * E.mean()).tofile(tmp.name)
import shutil; shutil.copy(real[:-4] + ".json", tmp.name[:-4] + ".json")   # the sidecar travels with the field
run = lambda **env: subprocess.run([sys.executable, os.path.join(HERE, "study2_run.py")], capture_output=True,
                                   text=True, env=dict(os.environ, **env)).stdout
outs = [run(STUDY2_REF_BIN=tmp.name) for _ in range(2)]
ok = all("=== VERDICT: REFUSED-GATED ===" in o and "GATE G7b dose-field hash: FAIL" in o for o in outs)
fails += not ok
print(f"  {'refused ' if ok else 'ACCEPTED'} garbage field, no control flag, two runs -> G7b both times")
# G7b, second attack: no pin file at all. The gate must refuse, not pin the field it just read.
o = run(STUDY2_REF_BIN=tmp.name, STUDY2_PINS=tmp.name[:-4] + ".pins-absent.json")
ok = "=== VERDICT: REFUSED-GATED ===" in o and "GATE G7b dose-field hash: FAIL" in o and "NONE" in o
fails += not ok
print(f"  {'refused ' if ok else 'ACCEPTED'} garbage field with no pin file -> G7b refuses instead of self-pinning")
for f in (tmp.name, tmp.name[:-4] + ".json"): os.unlink(f)
if sha() != before or sha(PINS) != pins_before:
    fails += 1
    print("  ACCEPTED a gated run rewrote verdict.json or pins.json")
else:
    print("  refused  gated runs left verdict.json and pins.json byte-identical")
print("all controls refused" if not fails else f"REFUSED: {fails} control(s) accepted")
sys.exit(1 if fails else 0)
