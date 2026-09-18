"""Render loop, ParaView 6.2 API: camera orbit + timestep toggle.

Cues are created via GetAnimationScene().Cues / servermanager; keyframes are
set through KeyFrames lists built from CameraKeyFrame with KeyTime/Position/
FocalPoint/ViewUp (the only properties this build exposes).
"""
import math
import os
import paraview.simple as smp

OUT = os.path.dirname(os.path.abspath(__file__))     # paraview/
FRAMES = f"{OUT}/frames"
os.makedirs(FRAMES, exist_ok=True)
for f in os.listdir(FRAMES):
    if f.endswith(".png"):
        os.remove(os.path.join(FRAMES, f))

smp._DisableFirstRenderCamera = True

pvd = smp.PVDReader(FileName=f"{OUT}/prrt_cycles.pvd")
v = smp.GetActiveViewOrCreate("RenderView")
v.ViewSize = [1600, 900]
v.Background = [0.04, 0.05, 0.08]
try:
    v.OrientationAxesVisibility = 0
except Exception:
    pass

disp = smp.Show(pvd, v)
disp.Representation = "Surface"
smp.ColorBy(disp, ("POINTS", "dose_Gy"))
lut = smp.GetColorTransferFunction("dose_Gy")
lut.RescaleTransferFunction(0.0, 25.0)
for preset in ("Inferno", "inferno", "Plasma", "Rainbow Uniform", "Cool to Warm"):
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
scene.NumberOfFrames = 120

N_ORBIT = 120
orbit = smp.KeyFrameAnimationCue()
orbit.AnimatedPropertyName = "Center"
orbit.AnimatedProxy = v
c_kf = []
for i in range(N_ORBIT):
    th = 2.0 * math.pi * i / (N_ORBIT - 1)
    r = 1.9
    kf = smp.CameraKeyFrame()
    kf.KeyTime = i / (N_ORBIT - 1)
    kf.Position = [0.5 + r * math.cos(th), 0.5 + r * math.sin(th), 1.35]
    kf.FocalPoint = [0.5, 0.5, 0.5]
    kf.ViewUp = [0, 0, 1]
    c_kf.append(kf)
orbit.KeyFrames = c_kf
scene.Cues = list(scene.Cues) + [orbit]

# timestep toggle via a proxy-property cue on the reader's TimestepValues
ts_cue = smp.KeyFrameAnimationCue()
ts_cue.AnimatedPropertyName = "TimestepValues"
ts_cue.AnimatedProxy = pvd
kf0 = smp.KeyFrame()
kf0.KeyTime = 0.0
kf0.KeyValues = [0]
kf1 = smp.KeyFrame()
kf1.KeyTime = 0.5
kf1.KeyValues = [0]
kf2 = smp.KeyFrame()
kf2.KeyTime = 0.5001
kf2.KeyValues = [1]
kf3 = smp.KeyFrame()
kf3.KeyTime = 1.0
kf3.KeyValues = [1]
ts_cue.KeyFrames = [kf0, kf1, kf2, kf3]
scene.Cues = list(scene.Cues) + [ts_cue]

smp.Render()
smp.SaveAnimation(f"{FRAMES}/frame.png", v, ImageResolution=[1600, 900],
                  FrameRate=24)
n = len([f for f in os.listdir(FRAMES) if f.endswith(".png")])
print("frames written:", n)