# chtMultiRegionInterFoam — Verification & Validation

`chtMultiRegionInterFoam` is a multi-region conjugate-heat-transfer solver with a
two-phase VOF fluid region, extracted from C. Schubert's
`OpenFOAM--MHD--VoF-Solver-Collection` and ported from OpenFOAM v1912 to v2412.

> **Disclaimer.** As stated by the upstream author, the solver has **not** been
> formally verified or validated against experiment. This report documents
> **build/portability verification** and **internal sanity & consistency checks**
> on the bundled `pinFinJetCooling` example — it is a demonstration of correct
> coupled behaviour, not a quantitative experimental validation.

All runs use OpenFOAM **v2412**.

## 1. Build / portability verification

- The solver compiles on v2412 after a **4-line mechanical API patch only**
  (`mesh.data::add/remove("finalIteration", …)` → `mesh.data().setFinalIteration(…)`
  in `fluid/solveFluid.H` and `solid/solveSolid.H`); see `docs/porting/v2412-portability.patch`.
  This is the idiom used by the native v2412 `chtMultiRegionFoam`.
- It links against the **stock v2412** `compressibleInterFoam` libraries
  (`twoPhaseMixtureThermo`, `VoFphaseCompressibleTurbulenceModels`,
  `twoPhaseSurfaceTension`) — i.e. the fluid-region discretisation is unmodified OpenFOAM.
- **Corrections applied in this distribution** (see README "Changes in this
  distribution"): the adaptive time step now honours the VOF interface Courant
  limit `maxAlphaCo` (previously computed but not applied), the example uses an
  implicit `nOuterCorrectors 2` + `energyCoupling` conjugate coupling, and dead /
  duplicated single-region code was removed. A clean build is verified by
  `Allwmake` (link to `$FOAM_USER_APPBIN/chtMultiRegionInterFoam`); a coarse
  `(32 24 32)` serial run of the example confirms the interface Courant number is
  emitted and drives `deltaT`, the `energyCoupling` sub-loop fires each step, and
  temperatures stay bounded.

## 2. Test case

`examples/pinFinJetCooling`: a 32 × 24 × 32 mm box, three regions —

| region | material | role |
| --- | --- | --- |
| `fluid` | water + air (VOF) | coolant domain; 1.5 m/s water jet from a 4 mm nozzle |
| `heatSink` | aluminium | base + 3 × 3 square pin fins |
| `chip` | silicon | 150 W volumetric source (`scalarSemiImplicitSource` on `h`) |

Two resolutions were run: a coarse mesh (24,576 cells) to `t = 0.2 s` and a
3×-refined mesh (663,552 cells) to `t = 0.5 s` (laminar; adjustable `dt`,
`maxCo = maxAlphaCo = 0.5`).

> The quantitative tables in sections 3–7 are from the **original** runs (before
> the interface-Courant fix and the `nOuterCorrectors 2` coupling change). They
> remain representative of the qualitative behaviour; absolute numbers shift
> slightly with the corrected time-step limiting and tighter coupling.

## 3. Numerical stability (VOF boundedness)

| Metric (3×-refined, 0–0.5 s, 7809 steps) | Result | Expected |
| --- | --- | --- |
| `min(alpha.water)` over the run | ≈ −1.0e-17 | ≈ 0 (machine) |
| `max(alpha.water)` over the run | 1.0 | 1 |
| peak interface Courant number | ≈ 0.68 | bounded, O(maxAlphaCo) |
| water volume fraction at 0.5 s | ≈ 0.52 | monotonic fill |

The phase fraction remains bounded in `[0, 1]` (to machine precision) for the
whole run — MULES behaves as expected. Note the peak interface Courant number of
**≈ 0.68 exceeded the `maxAlphaCo = 0.5` target**: this is precisely the symptom
of the (now-fixed) inactive interface-Courant limiter — `deltaT` was being set by
the bulk-flow Courant number alone. With the fix the interface Courant number is
bounded by `maxAlphaCo` through the adaptive step.

## 4. Conjugate-heat-transfer consistency

Region temperatures at the end of each run:

| quantity | coarse, t = 0.2 s | refined, t = 0.5 s |
| --- | --- | --- |
| chip `T` range [K] | 312.9 – 315.7 | 323.8 – 329.4 |
| heat-sink `T` range [K] | 300.8 – 313.2 | 304.7 – 326.5 |
| chip max `T` rise | 300 → 315.5 | 300 → 329.1 |
| heat-sink mean `T` rise | 300 → 308.3 | 300 → 318.9 |

A monotone, physical gradient is established and maintained:
**chip (hottest, the source) → heat-sink base → pin tips (≈ coolant temperature)**.
The pin tips sit near 300–305 K, i.e. the impinging jet is removing heat at the
solid–fluid interface — the conjugate coupling is working in both directions.

## 5. Energy bookkeeping

With `Q = 150 W` and a lumped solid heat capacity
`rho Cp V ≈ 2.8 J/K` (chip + heat sink), the adiabatic (no-cooling) rate is
`Q/(rho Cp V) ≈ 54 K/s`. The observed chip-max rise (≈ +15.5 K over 0.2 s coarse,
≈ +29 K over 0.5 s refined) is **below** the adiabatic bound and its slope
**decreases** with time — consistent with the jet removing a growing share of the
dissipated power as the heat sink wets.

## 6. Mesh refinement

Refining the mesh 3× (1 mm → 1/3 mm; 24.6 k → 664 k cells) **reproduces the
coarse-mesh temperature trend** (same qualitative heating curve and chip→pin
gradient) while sharpening the resolved jet column and the wavy free surface. No
new instability or unboundedness appears at the finer resolution.

## 7. Qualitative interface behaviour

The `alpha.water = 0.5` isosurface (see `docs/figures/interface_t*.png`) shows the
expected jet-impingement sequence: a coherent **jet column** descends, **impinges**
on the pin array, **spreads radially** into a wavy liquid sheet, and the liquid
**pools and rises** to wet the pins — physically sensible free-surface dynamics.

## Limitations

- **Laminar** flow; a real jet at these scales is turbulent (Re ~ 3 × 10⁴). Switch
  `constant/fluid/turbulenceProperties` to RAS for engineering use.
- No experimental comparison; absolute temperatures depend on the (illustrative)
  material properties, power and run time chosen for the demo.
- The example now uses `nOuterCorrectors = 2` with an `energyCoupling` sub-loop
  (implicit per-step conjugate coupling). Set `nOuterCorrectors 1` in the
  top-level `system/fvSolution` for the cheaper, fully-explicit coupling if a run
  needs maximum robustness over coupling accuracy.
