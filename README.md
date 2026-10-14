# chtMultiRegionInterFoam

![OpenFOAM v2412](https://img.shields.io/badge/OpenFOAM-v2412-blue)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green)

A **multi-region conjugate-heat-transfer (CHT) solver with a two-phase VOF fluid region** for OpenFOAM. It combines, in a single executable, what `interFoam`/`compressibleInterFoam` and `chtMultiRegionFoam` do separately:

- **fluid region(s)** — two compressible, non-isothermal, immiscible phases (e.g. water + air) captured with an algebraic-MULES **volume-of-fluid** interface, with surface tension, buoyancy and an energy equation;
- **solid region(s)** — heat conduction, optionally with volumetric heat sources (`fvOptions`);
- coupled at every fluid–solid and solid–solid interface by a **conjugate (mixed temperature/heat-flux) boundary condition**, so a free-surface flow and the heat conduction in the surrounding solids are solved together.

The bundled example is **liquid-jet-impingement cooling of a pin-fin heat sink on a hot GPU chip**: a water jet falls through air onto an aluminium pin-fin heat sink mounted on a heat-dissipating silicon die.

| | |
|:---:|:---:|
| ![jet impinging on the pin-fin heat sink (t = 0.3 s)](docs/figures/interface_t0.3.png) | ![chip / heat-sink temperature history](docs/figures/temperature_history.png) |

> ⚠️ **Provenance & scope.** This is **not a newly derived solver**. It is an extraction of `chtMultiRegionInterFoam` from C. Schubert's GPL-3 [`OpenFOAM--MHD--VoF-Solver-Collection`](https://github.com/c-schubert/OpenFOAM--MHD--VoF-Solver-Collection) (developed for OpenFOAM **v1912**), ported to **v2412** with a **minimal, documented 4-line API patch** (`mesh.data::add/remove(...)` → `mesh.data().setFinalIteration(...)`, see [`docs/porting/`](docs/porting)). The core fluid/solid discretisation is faithful to upstream; this distribution additionally **fixes the VOF interface-Courant time-step limiting** and **removes dead/duplicated code** (see [Changes in this distribution](#changes-in-this-distribution)). As the upstream author notes, the solver has not been formally verified or validated against experiment — the example here is a **demonstration** with internal sanity/consistency checks (see [V&V](docs/verification_validation.md)). Use accordingly.

## Physics

Per time step the multi-region PIMPLE driver solves each region in turn and couples them through the interface temperature condition. With `nOuterCorrectors > 1` (the example default), the first outer corrector is followed by an **energy-coupling sub-loop** that re-solves the fluid temperature and solid conduction with the flow frozen, making the conjugate coupling implicit within the step; `nOuterCorrectors 1` falls back to a single, fully-explicit (segregated) fluid-then-solid update.

**Fluid region** (`compressibleInterFoam` formulation, phase 1 = liquid):

| Equation | Field |
| --- | --- |
| VOF phase fraction (MULES): `d(alpha1)/dt + div(alpha1 U) + div(alpha1 alpha2 Ur) = 0` | `alpha.water` |
| Mixture momentum (PIMPLE) with gravity + surface tension `sigma kappa grad(alpha1)` | `U`, `p_rgh`, `p` |
| Energy (temperature form), mixture `rho = alpha1 rho1 + alpha2 rho2` | `T` |
| Compressible turbulence (here `laminar`; RAS/LES selectable) | — |

**Solid region(s):**

| Equation | Field |
| --- | --- |
| Heat conduction `d(rho h)/dt - laplacian(alpha_h, h) = q_source` | `h` (→ `T`) |

**Interface coupling:** `compressible::turbulentTemperatureRadCoupledMixed` enforces temperature and heat-flux continuity (`T_s = T_f`, `k_s dT_s/dn = k_f dT_f/dn`) on the auto-generated `<region>_to_<neighbour>` patches.

## Changes in this distribution

Relative to the upstream v1912 source (beyond the mechanical v2412 portability patch), this distribution applies the following corrections and clean-ups. The fluid/solid finite-volume discretisation is otherwise unchanged.

- **Interface-Courant time-step limiting (bug fix).** The adaptive time-step controller (`include/setMultiRegionDeltaT.H`) previously limited `deltaT` only by the bulk-flow Courant number and solid diffusion number; the VOF interface Courant number (`maxAlphaCo`) was computed and printed but **never applied**, so `maxAlphaCo` in `controlDict` had no effect. The controller now limits the fluid step by `min(maxCo/CoNum, maxAlphaCo/alphaCoNum)`, matching `compressibleInterFoam`. The initial step (`setInitialMultiRegionDeltaT.H`) is limited the same way.
- **Conjugate coupling demonstrated.** The example now runs with `nOuterCorrectors 2` and an `energyCoupling` sub-loop (implicit per-step fluid<->solid coupling) instead of the fully-explicit `nOuterCorrectors 1`. Both modes are supported; see `system/fvSolution`.
- **Adaptive outer-corrector reduction (M1, optional, off by default).** The driver can drop a time step to a single outer corrector (no `energyCoupling` sub-loop; that corrector is the final iteration, so `p_rgh`, `T` and `h` all use their Final solver settings - this needed a one-line fix in `fluid/pEqn.H`, which now selects `p_rghFinal` by the driver's `finalIter` instead of the configured outer-corrector count) once the relative per-step change of the flow and of the solid thermal field has stayed below `adaptiveOuterTol` for `adaptiveOuterSteps` consecutive steps (PIMPLE keys in the top-level `system/fvSolution`; defaults `0` and `5`). With `adaptiveOuterTol 0` the solver is bit-identical to the plain fixed-corrector driver. It pays off in quasi-steady / long thermal-design runs (in the tuning study at tolerance 0.5: -46% wall-clock on the baseline settings and a further -40% on the tuned preset, with at most 0.015 K chip-T drift); the consecutive-steps hysteresis keeps it dormant through transients at a conservative tolerance (at the loose 0.5 used in the study it also fires during start-up). See the [tuning study](docs/tuning_study.md) before enabling it.
- **Dead code removed.** Unused single-region leftovers (`VoF/setDeltaT.H`, `VoF/setRDeltaT.H`, `VoF/alphaCourantNo.H`, `VoF/alphaEqnSubCycle.H`, `VoF/createAlphaFluxes.H`, `fluid/alphaEqn_bu.H`, `fluid/setMultiRegionRDeltaT.H`), the unused `initialMass` bookkeeping and the dead `nEnergyCorrectors` control were removed. The dead outer-corrector `p_rgh.storePrevIter()` loop (and its `fluid/storeOldFluidFields.H` header) was removed as well (QW4; `p_rgh` is never under-relaxed across outer correctors, so results are unchanged).
- **Source hygiene.** The application identity/description was corrected (it self-identified as `chtMultiRegionTwoPhaseEulerFoam`), and all non-ASCII / work-in-progress comments were converted to ASCII English.

> **Pressure reference.** Like its parent `compressibleInterFoam`, the fluid solver sets the pressure level from the compressible thermodynamics and the boundary conditions; it does **not** carry a closed-domain `pRefCell`/`pRefValue` reference. Each fluid region therefore needs a pressure-referencing boundary (e.g. the `totalPressure` outlet in the example) or compressible pinning.

## Repository structure

```text
chtMultiRegionInterFoam/
├── Allwmake / Allwclean                  build scripts (wmake src/chtMultiRegionInterFoam)
├── src/
│   └── chtMultiRegionInterFoam/          solver source (built by Allwmake)
│       ├── chtMultiRegionInterFoam.C     application + shared multi-region setup headers
│       ├── createFields.H  createMeshes.H  …
│       ├── Make/                         files / options
│       ├── fluid/                        fluid-region (compressibleInterFoam VOF) headers
│       ├── solid/                        solid-region conduction headers
│       ├── include/                      multi-region time-control headers
│       └── VoF/                          header-only algebraic-MULES VoF (alpha equation)
├── examples/
│   └── pinFinJetCooling/                 pin-fin heat sink, liquid jet impingement cooling
└── docs/
    ├── verification_validation.md        sanity / consistency report
    ├── tuning_study.md                   performance tuning study (tuned preset, M1)
    ├── porting/                          upstream provenance + v2412 portability patch
    └── figures/                          result figures
```

> **Note on the layout.** `chtMultiRegionInterFoam` is built on a shared multi-region framework plus a bundled header-only VoF library. These are collected under a single solver directory `src/chtMultiRegionInterFoam/`, mirroring the layout of the sibling solvers (`dacFoam`, `adsorptionFoam`, `electrostaticSimpleFoam`); the multi-region sources are split into `fluid/`, `solid/` and `include/` subdirectories following the OpenFOAM `chtMultiRegionFoam` idiom. `Allwmake` builds this one application; the bundled `VoF/`, `fluid/`, `solid/` and `include/` headers are pulled in via the relative include paths in `Make/options`.

## Quick start

### Prerequisites

- OpenFOAM **v2412** (ESI/OpenCFD). The solver targets the bundled `compressibleInterFoam` libraries (`twoPhaseMixtureThermo`, `VoFphaseCompressibleTurbulenceModels`, `surfaceTensionModels`), which ship with v2412.
- For the result figures: ParaView (`pvbatch`); optional `gnuplot` for the temperature plot.

### Build

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc   # adjust to your installation
git clone https://github.com/swtbkim/chtMultiRegionInterFoam.git
cd chtMultiRegionInterFoam
./Allwmake
```

`./Allwmake` builds `chtMultiRegionInterFoam` into `$FOAM_USER_APPBIN`. Verify:

```bash
chtMultiRegionInterFoam -help
```

### Run the example

```bash
cd examples/pinFinJetCooling
./Allrun            # blockMesh -> topoSet -> createPatch -> splitMeshRegions -> decomposePar -> (parallel) solve -> reconstruct
```

The example is a 3-region case (`fluid`, `heatSink`, `chip`) and runs **in parallel** (16 ranks by default; see `system/decomposeParDict`). For a quick single-core run, coarsen the mesh in `system/blockMeshDict` (e.g. `(32 24 32)`) and replace the parallel block of `Allrun` with a plain `runApplication $(getApplication)`. Clean with `./Allclean`. See [`examples/pinFinJetCooling/README.md`](examples/pinFinJetCooling/README.md).

## Verification & validation

This is a demonstration solver; the checks below are internal sanity/consistency tests, **not** experimental validation (full report: [docs/verification_validation.md](docs/verification_validation.md)):

- **Builds** on OpenFOAM v2412 with only the documented 4-line portability patch; links against the stock `compressibleInterFoam` libraries.
- **Bounded VOF**: across the full 0.5 s run the phase fraction stays bounded (`min(alpha.water) ≈ -1e-17`, `max = 1`); the interface Courant number is now actively limited by `maxAlphaCo` through the adaptive time step.
- **Conjugate heat transfer is consistent**: with a 150 W chip source, the chip heats from 300 K toward a quasi-steady level while the jet holds the pin tips near the 300 K coolant temperature, giving a physical chip → base → pin-tip gradient.
- **Energy bookkeeping**: the chip-temperature rise rate is consistent with the lumped `Q / (rho Cp V)` estimate minus the (growing) jet heat removal.
- **Mesh refinement**: a 3×-refined mesh (664 k cells) reproduces the coarse-mesh temperature trend and sharpens the resolved jet/interface.

## Documentation

- [Verification & validation report](docs/verification_validation.md)
- [Performance tuning study](docs/tuning_study.md) — baseline vs tuned settings with per-step measurements, the opt-in tuned preset (-28% wall-clock, accuracy-neutral) and the M1 adaptive outer-corrector validation.
- [Step-by-step tutorial](https://github.com/swtbkim/openfoam-tutorial-ko) (Korean, 예제 19 of `openfoam-tutorial-ko`) — physics, file-by-file case setup, meshing, running (serial & parallel) and post-processing, with the jet-interface and temperature-history results, plus an advanced-options / performance-tuning chapter that walks through each tuning-tested setting and the M1 option.

## Citation

If you use this solver, please cite it via [CITATION.cff](CITATION.cff) **and** acknowledge the upstream collection (below).

## Acknowledgements and license

`chtMultiRegionInterFoam` is derived from **C. Schubert, [`OpenFOAM--MHD--VoF-Solver-Collection`](https://github.com/c-schubert/OpenFOAM--MHD--VoF-Solver-Collection)** (GPL-3), which itself combines the OpenFOAM standard solvers `compressibleInterFoam` and `chtMultiRegionFoam` (Copyright OpenFOAM Foundation / OpenCFD Ltd) with ideas from the DLR VoFLibrary, `varRhoTurbVOF` (Fan & Anglart) and isoAdvector (Roenby et al.). Please credit those works. The original upstream README is preserved in [`docs/porting/UPSTREAM_README.md`](docs/porting/UPSTREAM_README.md).

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via [openfoam.com](https://openfoam.com), and owner of the OPENFOAM and OpenCFD trade marks.

Distributed under the [GNU General Public License v3](LICENSE), as required for OpenFOAM derivatives.
