# pinFinJetCooling — liquid jet-impingement cooling of a pin-fin heat sink on a GPU chip

A three-region conjugate-heat-transfer + VOF example: a downward **water jet**
falls through air onto an aluminium **pin-fin heat sink** (base + 3 × 3 pins)
mounted on a heat-dissipating silicon **chip**.

| region | material | notes |
| --- | --- | --- |
| `fluid` | water + air (VOF) | 1.5 m/s jet from a 4 mm central `nozzle`, atmospheric `outlet` |
| `heatSink` | aluminium (k = 205) | base + 3 × 3 square pins |
| `chip` | silicon (k = 149) | 150 W volumetric source (`system/chip/fvOptions`) |

The single background mesh is split into the three region meshes with
`splitMeshRegions`, which auto-creates the coupled `<region>_to_<neighbour>`
interfaces. Initial/boundary fields are shipped per region in `0.orig/`.

## Run

```bash
./Allrun       # blockMesh -> topoSet -> createPatch -> splitMeshRegions
               # -> decomposePar(-allRegions) -> parallel solve -> reconstructPar
./Allclean     # reset to a clean case
```

`Allrun` runs **in parallel** (4 ranks; `system/decomposeParDict` and the
per-region `system/<region>/decomposeParDict`).

For a **quick serial run**, coarsen `system/blockMeshDict` to `(32 24 32)` and
replace the decompose/parallel/reconstruct lines of `Allrun` with:

```bash
runApplication $(getApplication)
```

## Region coupling

The top-level `system/fvSolution` sets `nOuterCorrectors 1`: the cheapest,
fully-explicit coupling (one fluid solve then one solid solve per step). For
fluid<->solid conjugate coupling solved **implicitly within each time step**,
set `nOuterCorrectors 2` there and add an `energyCoupling { iterations 2; }`
sub-dictionary to `PIMPLE`. The VOF interface
Courant limit `maxAlphaCo` (in `system/controlDict`) is honoured by the adaptive
time step.

## Expected results

(Indicative; from the original `nOuterCorrectors 1` refined run.)

- The water jet impinges on the pins, spreads, and pools (the `alpha.water = 0.5`
  interface; ~52 % fill at 0.5 s on the refined mesh).
- The chip heats from 300 K (≈ 329 K at 0.5 s) while the jet holds the pin tips
  near the 300 K coolant — a chip → base → pin-tip temperature gradient.
- Monitors (printed each step / written under `postProcessing/`):
  `chipTemperature` (max `T` in the chip) and `heatSinkTemperature` (mean `T`).

## Post-processing

```bash
paraFoam -touch-all       # then open in ParaView (Contour alpha.water = 0.5; T colour map)
# or reproduce the report figures (interface isosurface frames):
foamToVTK -allRegions -time '0,0.1,0.2,0.3,0.4,0.5' -fields '(alpha.water T)'
pvbatch --force-offscreen-rendering render_interface.py   # writes figures/interface_t<time>.png
```

A step-by-step walkthrough of this case (Korean) is example 19 of
[openfoam-tutorial-ko](https://github.com/swtbkim/openfoam-tutorial-ko).
