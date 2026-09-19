#!/usr/bin/env python3
"""Plant defects in a scratch copy of the tree and run the untouched runner there.

Each case must end REFUSED-GATED at the gate named here, or the gate is
decoration (issues #11, #13). The shipped runner carries no branch that exists
only for testing: no flag, no environment override, no control switch. Every
defect is a file mutation in the copy, and the mutation is checked to have
applied before the run. The first check below is that greppable criterion.

    python3 study2_response_layer/controls.py     # exit 0 or 1
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
RUNNER = os.path.join(HERE, "study2_run.py")
FILES = ["study2_response_layer/study2_run.py", "study2_response_layer/pins.json",
         "study2_response_layer/verdict.json", "study1/runs/ref_p40_h3.bin",
         "study1/runs/ref_p40_h3.json", "study1/geometry_pitch40_meta.json",
         "data/geometry_pitch40.npz"]
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
fails = 0
def report(ok, text):
    global fails; fails += not ok; print(f"  {'refused ' if ok else 'ACCEPTED'} {text}")

# 0. the stopping criterion: no test-only branch in the shipped runner
src = open(RUNNER, encoding="utf-8").read()
report(not re.search(r"STUDY2_|--control|--pin|\bCONTROL\b", src),
       "shipped runner has no test-only flag, override or switch")

def fresh():
    d = tempfile.mkdtemp(prefix="study2-control-")
    for f in FILES:
        os.makedirs(os.path.dirname(os.path.join(d, f)), exist_ok=True)
        shutil.copy(os.path.join(ROOT, f), os.path.join(d, f))
    return d
def run(d):
    return subprocess.run([sys.executable, os.path.join(d, "study2_response_layer", "study2_run.py")],
                          capture_output=True, text=True).stdout
def refused_at(out, gate):
    return "=== VERDICT: REFUSED-GATED ===" in out and f"GATE {gate}: FAIL" in out
def write_field(d, arr):
    p = os.path.join(d, "study1/runs/ref_p40_h3.bin"); before = sha(p); arr.tofile(p); assert sha(p) != before
    return sha(p)
def repin(d, h):   # a planted field with a matching pin, so a later gate is the one under test
    p = os.path.join(d, "study2_response_layer/pins.json"); j = json.load(open(p)); j["ref_p40_h3.bin"] = h
    json.dump(j, open(p, "w"), indent=1)

E = np.fromfile(os.path.join(ROOT, "study1/runs/ref_p40_h3.bin"), dtype=np.float64)
g = np.load(os.path.join(ROOT, "data/geometry_pitch40.npz")); viable = (g["cell_id"] > 0).ravel()

# 1. G6b: a flat field with 0.71% MC noise, pinned so only G6b can refuse it
d = fresh(); h = write_field(d, np.full(E.size, float(E[viable].mean())) * (1.0 + 0.0071 * np.random.default_rng(1).standard_normal(E.size)))
repin(d, h); report(refused_at(run(d), "G6b noise-only null"), "flat field, 0.71% noise, matching pin -> G6b"); shutil.rmtree(d)

# 2. G8: the verdict's own uniform arm built from the median, a source mutation
d = fresh(); p = os.path.join(d, "study2_response_layer/study2_run.py"); s = open(p).read()
old = "    D_u = float(D_h.mean())\n"; assert s.count(old) == 1; open(p, "w").write(s.replace(old, "    D_u = float(np.median(D_h))\n"))
report(refused_at(run(d), "G8 energy bookkeeping"), "uniform arm from the median inside arm_pair -> G8"); shutil.rmtree(d)

# 3. G7b: a garbage field, two runs; run 2 refuses only if run 1 wrote nothing
d = fresh(); write_field(d, np.random.default_rng(2).uniform(0.1, 10.0, E.size) * E.mean())
v, pn = [os.path.join(d, "study2_response_layer", f) for f in ("verdict.json", "pins.json")]; hv, hp = sha(v), sha(pn)
outs = [run(d), run(d)]
report(all(refused_at(o, "G7b dose-field hash") for o in outs), "garbage field, two runs -> G7b both times")
report(sha(v) == hv and sha(pn) == hp, "gated runs left verdict.json and pins.json byte-identical")

# 4. G7b: the same field with pins.json and verdict.json deleted; no self-pin
os.unlink(v); os.unlink(pn); o = run(d)
report(refused_at(o, "G7b dose-field hash") and "NONE" in o and not os.path.exists(v),
       "garbage field, no pin file, no verdict file -> G7b refuses and writes nothing"); shutil.rmtree(d)

print("all controls refused" if not fails else f"REFUSED: {fails} control(s) accepted")
sys.exit(1 if fails else 0)
