# Render the water-air interface for every time exported by foamToVTK.
# Usage (from the case directory, after foamToVTK -allRegions ...):
#   pvbatch --force-offscreen-rendering render_interface.py [caseDir]
import json, os, sys
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
NAME = os.path.basename(CASE)          # foamToVTK prefixes files with the case name
VTK  = CASE + "/VTK"
OUT  = CASE + "/figures"
NONE = ['POINTS', '']
os.makedirs(OUT, exist_ok=True)

def time_label(t):
    return "%.1f" % t if round(t, 1) == t else "%g" % t

# (vtm file name, time label) for each exported time, from the fluid series file
with open(f"{VTK}/fluid/{NAME}.vtm.series") as f:
    FRAMES = [(e["name"], time_label(e["time"])) for e in json.load(f)["files"]]

def internal(vtmfile):
    r  = XMLMultiBlockDataReader(FileName=[vtmfile])
    eb = ExtractBlock(Input=r); eb.Selectors = ['//internal']
    m  = MergeBlocks(Input=eb)
    return r, eb, m

def render_frame(vtm, tlabel):
    fr, feb, fluidIM = internal(f"{VTK}/fluid/{vtm}")
    hr, heb, hsIM    = internal(f"{VTK}/heatSink/{vtm}")
    cr, ceb, chipIM  = internal(f"{VTK}/chip/{vtm}")

    # water-air interface alpha.water = 0.5
    fluidP = CellDatatoPointData(Input=fluidIM)
    interface = Contour(Input=fluidP, ContourBy=['POINTS', 'alpha.water'], Isosurfaces=[0.5])

    # domain outline = feature edges of the clean (manifold) fluid surface
    fluidS = ExtractSurface(Input=fluidIM)
    edges  = FeatureEdges(Input=fluidS)
    edges.FeatureAngle = 30
    edges.FeatureEdges = 1
    edges.BoundaryEdges = 1
    edges.NonManifoldEdges = 0
    edges.ManifoldEdges = 0

    rv = CreateRenderView()
    rv.ViewSize = [1100, 850]
    rv.Background = [1, 1, 1]
    rv.OrientationAxesVisibility = 0
    rv.UseColorPaletteForBackground = 0

    dh = Show(hsIM, rv);   dh.Representation = 'Surface'; dh.ColorArrayName = NONE
    dh.DiffuseColor = [0.72, 0.72, 0.75]; dh.Ambient = 0.30
    dc = Show(chipIM, rv); dc.Representation = 'Surface'; dc.ColorArrayName = NONE
    dc.DiffuseColor = [0.90, 0.35, 0.10]; dc.Ambient = 0.30

    di = Show(interface, rv); di.ColorArrayName = NONE
    di.DiffuseColor = [0.10, 0.45, 0.95]; di.Ambient = 0.20

    de = Show(edges, rv); de.ColorArrayName = NONE
    de.DiffuseColor = [0, 0, 0]; de.LineWidth = 1.4

    txt = Text(); txt.Text = f"t = {tlabel} s"
    dt = Show(txt, rv); dt.Color = [0, 0, 0]; dt.FontSize = 22
    dt.WindowLocation = 'Upper Left Corner'

    rv.ResetCamera()
    cam = GetActiveCamera()
    cam.Azimuth(35); cam.Elevation(22)
    rv.ResetCamera()
    cam.Dolly(1.1); rv.Update()

    Render(rv)
    fn = f"{OUT}/interface_t{tlabel}.png"
    SaveScreenshot(fn, rv, ImageResolution=[1100, 850], TransparentBackground=0)
    print("saved", fn)

    for o in (dt, txt, de, edges, fluidS, di, interface, fluidP,
              dc, chipIM, ceb, cr, dh, hsIM, heb, hr, fluidIM, feb, fr, rv):
        Delete(o)

for vtm, t in FRAMES:
    render_frame(vtm, t)
print("DONE")
