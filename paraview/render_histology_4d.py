"""Render the 4D histology loop: 4 cycles x camera orbit, 2 fields."""
import math
import os
import paraview.simple as smp

OUT = "/home/aurascoper/Developer/PRRT-spatial-cpm/study3_closed_loop/paraview"
FRAMES = f"{OUT}/frames_4d"
os.makedirs(FRAMES, exist_ok=True)
for f in os.listdir(FRAMES):
    if f.endswith(".png"):
        os.remove(os.path.join(FRAMES, f))

smp._DisableFirstRenderCamera = True

pvd = smp.PVDReader(FileName=f"{OUT}/histology_4d.pvd")
v = smp.GetActiveViewOrCreate("RenderView")
v.ViewSize = [1600, 900]
v.Background = [0.04, 0.05, 0.08]
try:
    v.OrientationAxesVisibility = 0
except Exception:
    pass

disp = smp.Show(pvd, v)
disp.Representation = "Surface"
smp.ColorBy(disp, ("POINTS", "ln_expression"))
lut = smp.GetColorTransferFunction("ln_expression")
lut.RescaleTransferFunction(-1.5, 1.5)
for preset in ("Viridis", "viridis", "Plasma", "Cool to Warm"):
    try:
        lut.ApplyPreset(preset, True)
        break
    except RuntimeError:
        continue

v.ResetCamera()
v.CameraPosition = [1.9, 1.9, 1.9]
v.CameraFocalPoint = [0.5, 0.5, 0.5]
v.CameraViewUp = [0, 0, 1]

scene = smp.GetAnimationScene()
scene.PlayMode = "Sequence"
scene.NumberOfFrames = 160

# cut away a quarter so the interior histology is visible
clip = smp.Clip(Input=pvd)
clip.ClipType = "Box"
clip.Invert = 1
clip.ClipType.Bounds = [0.5, 20.0, 0.5, 20.0, 0.5, 20.0]
smp.Hide(pvd, v)
disp2 = smp.Show(clip, v)
smp.ColorBy(disp2, ("POINTS", "ln_expression"))
lut2 = smp.GetColorTransferFunction("ln_expression")
lut2.RescaleTransferFunction(-3.0, 1.0)
try:
    lut2.ApplyPreset("Viridis", True)
except RuntimeError:
    pass

# orbit cue
import math
N = 96
orbit = smp.KeyFrameAnimationCue()
orbit.AnimatedPropertyName = "Center"
orbit.AnimatedProxy = v
kfs = []
for i in range(N):
    th = 2.0 * math.pi * i / (N - 1)
    r = 1.7
    kf = smp.CameraKeyFrame()
    kf.KeyTime = i / (N - 1)
    kf.Position = [0.5 + r * math.cos(th), 0.5 + r * math.sin(th), 1.15]
    kf.FocalPoint = [0.5, 0.5, 0.5]
    kf.ViewUp = [0, 0, 1]
    kfs.append(kf)
orbit.KeyFrames = kfs
scene.Cues = list(scene.Cues) + [orbit]

# timestep cue: quarter of the loop per cycle
ts_cue = smp.KeyFrameAnimationCue()
ts_cue.AnimatedPropertyName = "TimestepValues"
ts_cue.AnimatedProxy = pvd
ts_kfs = []
for i, ts in enumerate([0, 0, 1, 1, 2, 2, 3, 3]):
    kf = smp.KeyFrame()
    kf.KeyTime = i / 7.0
    kf.KeyValues = [float(i // 2)]
    ts_kf = kf
    if i == 0:
        ts_kf.KeyTime = 0.0
    ts_kfs.append(ts_kf)
# build explicit keyframes: hold each cycle for 1/4 of the loop
ts_cue.KeyFrames = []
for idx, (t, val) in enumerate([(0.0, 0), (0.25, 0), (0.251, 1), (0.5, 1),
                                 (0.501, 2), (0.75, 2), (0.751, 3), (1.0, 3)]):
    kf = smp.KeyFrame()
    kf.KeyTime = t
    kf.KeyValues = [float(val)]
    ts_cue.KeyFrames.append(kf)
scene.Cues = list(scene.Cues) + [ts_cue]

smp.Render()
smp.SaveAnimation(f"{FRAMES}/frame.png", v, ImageResolution=[1600, 900],
                  FrameRate=24)
n = len([f for f in os.listdir(FRAMES) if f.endswith(".png")])
print("frames written:", n)