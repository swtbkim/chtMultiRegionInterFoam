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

`Allrun` runs **in parallel** (16 ranks; `system/decomposeParDict` and the
per-region `system/<region>/decomposeParDict`). On a host with fewer physical
cores than ranks, `Allrun` exports `OMPI_MCA_rmaps_base_oversubscribe=1` so MPI
still starts (equivalently use `mpirun --use-hwthread-cpus`).

For a **quick serial run**, coarsen `system/blockMeshDict` to `(32 24 32)` and
replace the decompose/parallel/reconstruct lines of `Allrun` with:

```bash
runApplication $(getApplication)
```

## Region coupling

The top-level `system/fvSolution` sets `nOuterCorrectors 2` with an
`energyCoupling` sub-loop, so the fluid<->solid conjugate coupling is solved
**implicitly within each time step**. For the cheapest (fully-explicit,
fluid-then-solid) coupling, set `nOuterCorrectors 1` there. The VOF interface
Courant limit `maxAlphaCo` (in `system/controlDict`) is honoured by the adaptive
time step.

## Baseline vs tuned settings

The shipped dictionaries are the **conservative baseline** (`nCorrectors 2`,
`nAlphaSubCycles 2`, `p_rgh relTol 0.01`, upwind `div(rhoPhi,T)`, `scotch`
fluid decomposition) with the adaptive outer-corrector option (M1,
`adaptiveOuterTol` / `adaptiveOuterSteps` in `system/fvSolution`) present but
off. A **tuned preset** validated in the tuning study (about -28% wall-clock, accuracy-neutral)
and the M1 option are described in
[`docs/tuning_study.md`](../../docs/tuning_study.md); apply them by hand when needed.

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
