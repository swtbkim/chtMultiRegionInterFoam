# chtMultiRegionInterFoam - performance tuning study (pinFinJetCooling)

This report documents a step-by-step tuning study of performance and accuracy
levers for `chtMultiRegionInterFoam` (VOF two-phase fluid + multi-region
conjugate heat transfer) on the bundled `pinFinJetCooling` example. Each
candidate change was tuning-tested against the unmodified solver and example on
an identical case and kept only where the measurements justified it.

The example in `examples/pinFinJetCooling` deliberately keeps its conservative
baseline settings (`nCorrectors 2`, `nAlphaSubCycles 2`, `p_rgh relTol 0.01`,
upwind energy advection, `scotch` decomposition, M1 off); the accepted case
settings are an **opt-in preset** ([Tuned preset](#tuned-preset-opt-in)). Two
**code** items from the study are in the solver (`src/chtMultiRegionInterFoam/`):

- **QW4** - removal of a dead `p_rgh.storePrevIter()` loop (and the
  `fluid/storeOldFluidFields.H` header it used). No behavioural change.
- **M1** - optional adaptive outer-corrector reduction
  ([M1](#m1---adaptive-outer-corrector-reduction-off-by-default)). **Off by
  default** (`adaptiveOuterTol 0`), in which case results are identical to the
  solver without it.

A Korean step-by-step walkthrough of this study is example 19 of the tutorial
repository [openfoam-tutorial-ko](https://github.com/swtbkim/openfoam-tutorial-ko).

## Method

- **Baseline:** the shipped settings (`nOuterCorrectors 2` + `energyCoupling 2`,
  `DIC` smoother, `relTol 0.01`, `nAlphaSubCycles 2`, upwind energy advection).
  Each candidate is applied on top of the previously accepted ones (M1 off
  unless stated).
- **Test:** identical mesh template (blockMesh `(48 36 48)`, ~83k cells,
  `fluid` ~95%), serial, fixed `endTime = 0.02 s` (181 baseline steps), minimal
  I/O (`writeControl timeStep; writeInterval 1000000`). Only the change under
  test differs; M6 (decomposition) was evaluated separately.
- **Metrics:** `ExecutionTime` (wall-clock), per-equation linear-solver
  iteration totals (cost proxy), `p_rgh` GAMG cycles/solve, chip and heat-sink
  max-T (the QoI; from the solid regions' `Min/max T` log line at the last
  step), and `alpha.water` bounds at the last step (interface integrity).
- **Gate:** a change is adopted only if it improves wall-clock (or accuracy) AND
  keeps the QoI within tolerance AND keeps `alpha.water` bounded.
- **Wall-clock noise:** ~0.5% within one session (baseline re-run 61.44 vs
  61.10 s), much larger across sessions (the same baseline measured 61.1 s,
  58.68 s and 54.12 s in three sessions). Compare wall-clock only within one
  session; iteration counts are deterministic.

> **Scale caveat.** The test is serial at 83k cells over a 0.02 s window that is
> almost entirely the start-up transient. It captures the per-step cost levers
> (pressure solver, correctors, sub-cycles) but NOT (a) production parallel
> scaling (663k cells / 16 ranks), where the GAMG and decomposition levers M5/M6
> pay off, nor (b) the quasi-steady behaviour of the adaptive coupling levers
> (QW3, M1). M1 at a conservative tolerance (0.1) stayed dormant in this window;
> the M1 speed-ups below come from the loose tolerance 0.5, which also fires
> during the start-up transient.

## Baseline

| steps | execTime | p_rgh solves / iters / cyc | T iters | h iters | chip Tmax | alpha |
|---|---|---|---|---|---|---|
| 181 | **61.1 s** | 724 / 2776 / 3.83 | 1230 | 898 | 302.985 K | -2.8e-11 .. 1 |

Per step the driver does 4 `p_rgh` solves (`nOuterCorrectors 2` x
`nCorrectors 2`), 4 fluid-T solves (1 + `energyCoupling` 2 + 1) and 8 solid-h
solves (4 per step x 2 regions); `U` is never solved (`momentumPredictor no`);
alpha is explicit MULES (0 linear iterations). The pressure solve dominates - it
is the optimization target.

## Per-step tuning-test results

Each step is compared with the configuration it was applied to. The
QW2 -> free -> QW5a -> QW3 chain gives the **intermediate configuration** on
which M5, M2, M3, QW5b and the M1 calibration runs were tested (the "54.3 s"
below): 183 steps, 732 `p_rgh` solves / 2336 iters (3.19 cycles/solve),
T 1951 iters, chip Tmax 302.9833 K, 54.28 s.

| Step | Change | Wall-clock | QoI (chip Tmax) | alpha | Verdict |
|---|---|---|---|---|---|
| QW1 | `p_rgh` smoother DIC -> DICGaussSeidel | 64.2 s vs 61.1 s: +5% (slower) | identical | ok | **reverted** (serial); production candidate |
| QW2 | intermediate `p_rgh` relTol 0.01 -> 0.05 | 58.6 s: -4% | +0.0002 K | ok | **adopted** |
| QW4 | delete dead `p_rgh.storePrevIter()` (code) | 0 (bit-identical) | identical | ok | **adopted**, in solver (hygiene) |
| free | `nAlphaSubCycles` 2 -> 1 | 52.4 s: -10% more | -0.007 K | ok | **adopted** |
| QW5a | `div(rhoPhi,T)` upwind -> limitedLinear 1 | 54.2 s: +3% (cost) | +0.005 K (2nd order) | ok | **adopted** (accuracy) |
| QW3 | `energyCoupling` `convergence {"T" 1e-3;}` | 54.3 s: neutral | identical | ok | **adopted** (can only skip the whole sub-loop; see below) |
| M1 | adaptive outer-corrector reduction (code) | off by default; tol 0.5: -46% on baseline, -40% on preset | -0.015 K / +0.0006 K | ok (preset + M1: -4.3e-8) | **in solver, off by default**; see [M1](#m1---adaptive-outer-corrector-reduction-off-by-default) |
| M5 | GAMG `nCellsInCoarsestLevel` 10 -> 200 | 55.3 s vs 54.3 s: +2% | -0.006 K | ok | **production-gated** (not in preset) |
| M6 | fluid decomposition scotch -> hierarchical | -3.5% proc faces | - | - | **in preset** (parallel only) |
| M2 | `nCorrectors` 2 -> 1 | 44.0 s vs 54.3 s: -19% | -0.003 K | ok | **in preset** (with caveat) |
| M3 | `maxAlphaCo` 0.5 -> 1.0 | 51.6 s vs 54.3 s: -5% | -0.008 K | overshoot 1.0000003 | **rejected** |
| QW5b | `div(rhoPhi,U)` upwind -> limitedLinearV 1 | 61.0 s vs 54.3 s: +12% (slower) | -0.003 K | overshoot 1.0000017 | **rejected** |
| M4 | cache frozen-branch contErr/divUp (code) | ~0 | not run | - | **rejected** (cheap explicit fvc ops, ~0 saving at 2 passes/step) |

### Notes on the instructive results

- **QW1 (reverted).** `DICGaussSeidel` cut GAMG cycles 5.5% (3.83 -> 3.62) but
  was ~5% slower: the heavier Gauss-Seidel sweep costs more than the cycle
  saving at 83k serial. Gate on wall-clock, not iterations. It remains a
  production candidate where the coarse-grid all-reduce cost dominates.
- **QW2 (adopted).** Loosening only the 3 intermediate `p_rgh` solves (keeping
  `p_rghFinal` tight at `relTol 0`) cut `p_rgh` iters 19% (2776 -> 2251) for -4%
  wall-clock at +0.0002 K.
- **QW3 (adopted, neutral here).** Same iteration counts as before it (366
  energy-coupling passes in 183 steps). In v2412 `loopControl::loop()` checks
  `convergence { "T" 1e-3; }` before every pass, including the first, against
  the initial residual of the step's **first** fluid-T solve (outer corrector 1,
  before the sub-loop). That value does not change between passes, so the
  criterion either skips the whole sub-loop (T 4 -> 2 and h 8 -> 4 solves per
  step) or never fires; with `iterations 2` it cannot give exactly 1 pass.
- **QW5a vs QW5b.** `div(rhoPhi,T)` (T is solved directly) is a clean 2nd-order
  accuracy upgrade at +3% cost. `div(rhoPhi,U)` is NOT cost-free: with
  `momentumPredictor no` its matrix still feeds `UEqn.A()/H()` -> `HbyA` -> the
  pressure equation, so it changed the pressure system, needed more steps
  (195 vs 183), ran +12% slower and grew the alpha overshoot. Rejected.
- **M2 (in preset, with caveat).** `nCorrectors` 2 -> 1 halved the pressure
  solves (732 -> 370) for -19% wall-clock at -0.003 K with alpha bounded. The
  naive "-50%" did NOT materialise: the single corrector starts from a worse
  guess, so cycles/solve rose 3.19 -> 5.14. Caveat: the 2nd corrector also
  stabilises the explicit surface-tension term; verify interface quality on a
  full-length run and restore `nCorrectors 2` if parasitic currents or alpha
  overshoot appear. This is one reason the shipped example keeps
  `nCorrectors 2`.
- **M3 (rejected).** After the interface-Courant fix (see README, "Changes in
  this distribution") `maxAlphaCo` is the binding time-step limiter (interface
  Co maxed at 0.499 vs bulk 0.417), so raising it to 1.0 did enlarge dt (165 vs
  183 steps). But the per-step cost rose, netting only -5%, and with a single
  sub-cycle the interface Courant number near 1.0 pushed alpha to 1.0000003
  (min -4.1e-7).

## Cumulative result (tuned preset, M1 off)

| | baseline | tuned preset | change |
|---|---|---|---|
| steps | 181 | 185 | +2% |
| execTime | 61.1 s | **44.11 s** | **-27.8%** |
| p_rgh solves | 724 | 370 | -49% |
| p_rgh iters | 2776 | 1902 | -31% |
| chip Tmax | 302.9849 K | 302.9807 K | -0.004 K (neutral) |
| heat-sink Tmax [1] | 301.4706 K | 301.4680 K | -0.003 K (neutral) |
| alpha.water | -2.8e-11 .. 1 | -1.4e-11 .. 1 | bounded |

[1] Region maximum from the heat-sink `Min/max T` log line at the last step.

With the tuned preset the solver runs **~1.39x faster** than with the shipped
baseline settings, accuracy-neutral (chip-T within 0.004 K), interface bounded,
and with a 2nd-order energy advection scheme. The win comes almost entirely
from halving the pressure-solve work (QW2 + M2) plus the cheaper alpha
sub-cycle.

## M1 - adaptive outer-corrector reduction (off by default)

After every step M1 measures the relative change of the flow
(`max|U - U.oldTime| / max|U.oldTime|` over the fluid regions) and of the solid
enthalpy `h` (same measure). When both stay below `adaptiveOuterTol` for
`adaptiveOuterSteps` **consecutive** steps, the next step runs a single outer
corrector instead of `nOuterCorrectors`, skipping the redundant pressure pass
and the energy-coupling sub-loop. The first step that fails the test resets the
counter, and the next step is fully coupled again.

```
PIMPLE                    // top-level system/fvSolution
{
    nOuterCorrectors    2;
    adaptiveOuterTol    0;    // 0 (default) = off, identical to the plain driver
    adaptiveOuterSteps  5;    // hysteresis: consecutive settled steps required
}
```

That single corrector is the final iteration (`finalIter` in the driver), so a
reduced step uses all Final solver settings (`p_rghFinal`, `TFinal`, `hFinal`)
and ends at the same final tolerances as a fully coupled step. To this end
`fluid/pEqn.H` selects `p_rghFinal` with the driver's `finalIter` instead of
the upstream `oCorr == nOuterCorr-1`; the two are equal whenever M1 does not
reduce a step. An earlier pre-release implementation used the upstream
expression, so reduced steps ended `p_rgh` at the intermediate tolerance. This
was fixed before release by selecting `p_rghFinal` with `finalIter`, which
costs pressure iterations (+638, +50%, on baseline + M1) but guarantees that
every reduced step ends at the Final tolerance.

When enabled, every step logs
`adaptiveOuter: flowChange=... thermChange=... steadyCount=... nextCorrectors=...`
and each reduced step starts with
`OUTER-EXIT: quasi-steady, using 1 outer corrector (tight)`. On a reduced step
the last `p_rgh` solve before the next `adaptiveOuter:` line (the second of two
with `nCorrectors 2`) is the `p_rghFinal` solve and must end below its absolute
tolerance (1e-7).

### Results

All four runs were made sequentially in one session with the current solver
(serial, 83k-cell template, 0.02 s), so their wall-clock times are comparable.
Percentages are relative to the first row; preset + M1 is -40.1% relative to
the preset alone.

| config | steps | outer-exits | execTime | p_rgh solves / iters / cyc | T / h iters | chip Tmax | heat-sink Tmax | alpha |
|---|---|---|---|---|---|---|---|---|
| baseline, M1 off | 181 | 0 | 54.12 s | 724 / 2776 / 3.83 | 1230 / 898 | 302.9849 | 301.47064 | -2.8e-11 .. 1 |
| baseline + M1 (tol 0.5, steps 5) | 177 | 156 | **29.38 s (-45.7%)** | 396 / 1910 / 4.82 | 534 / 802 | 302.97001 | 301.44507 | -7.3e-18 .. 1 |
| tuned preset, M1 off | 185 | 0 | 38.83 s (-28.3%) | 370 / 1902 / 5.14 | 1953 / 1006 | 302.9807 | 301.468 | -1.4e-11 .. 1 |
| tuned preset + M1 (tol 0.5, steps 5) | 185 | 150 | **23.27 s (-57.0%)** | 220 / 1775 / 8.07 | 891 / 859 | 302.98126 | 301.45326 | **-4.3e-8** .. 1 |

- **Baseline + M1:** 156 of 177 steps reduced; `p_rgh` solves 724 -> 396,
  fluid-T solves 724 -> 240, solid-h solves 1448 -> 480, energy-coupling passes
  362 -> 42. QoI drift vs M1 off: chip -0.015 K, heat-sink -0.026 K. Every
  reduced step's last `p_rgh` solve ended below 1e-7 (n = 156; min 3.35e-8,
  median 6.53e-8, max 9.99e-8).
- **Tuned preset + M1:** 150 of 185 steps reduced, -40.1% vs the preset; QoI
  drift vs the preset: chip +0.0006 K, heat-sink -0.015 K. Every reduced step's
  (only) `p_rgh` solve ended below 1e-7 (n = 150, max 9.98e-8). **Caveat:** the
  `alpha.water` undershoot grew to -4.3e-8 (vs ~1e-11 in the M1-off runs and
  ~1e-17 for baseline + M1): with `nCorrectors 1` a reduced step has a single
  pressure solve, which must reach `p_rghFinal` from a worse initial guess
  (8.07 cycles/solve). Still bounded in practice, but monitor `alpha.water` on
  full-length or surface-tension-dominated runs, and fall back to M1 off or
  `nCorrectors 2` if the undershoot grows.

### Calibration and hysteresis

M1's tolerance and hysteresis were calibrated in June with a **pre-release
build in which reduced steps still ended `p_rgh` at the intermediate
tolerance** (here QW2's `relTol 0.05`), on the intermediate configuration
(still `nCorrectors 2`). Use the [Results](#results) above for M1's current
cost and accuracy; this table justifies `adaptiveOuterSteps`. Speed-ups and QoI
drifts are relative to the M1-off row.

| config | steps | execTime | outer-exits | chip Tmax | alpha | note |
|---|---|---|---|---|---|---|
| off (`adaptiveOuterTol 0`) | 183 | 54.3 s | 0 | 302.9833 | -9.6e-10 .. 1 | reference; same iteration counts as without M1 |
| forced on (no warm-up) | 176 | 21.5 s | 174 | 302.9692 | -1.0e-17 .. 1 | -60%; chip -0.014 K, heat-sink -0.026 K |
| tol 0.5 + hysteresis (5 steps) | 180 | 28.8 s | 142 | 302.9795 | -1.5e-9 .. 1.0000001 | -47%; chip -0.004 K, heat-sink -0.018 K |
| tol 0.1, **no** hysteresis | 28 (crashed) | - | 3 | - | -3.3e-4 .. 1.0001 | **diverged** (fatal/nan in log) |
| tol 0.1 + hysteresis [2] | 183 | 54.8 s | 0 | 302.9833 | -9.6e-10 .. 1 | dormant; completes cleanly |

[2] This row's configuration is inferred from the study records; its exact
dictionary was not retained.

- **Hysteresis is required.** A naive tolerance is unsafe: an isolated
  low-change step during an active transient fired the reduction and diverged
  (alpha unbounded). The consecutive-steps hysteresis (`adaptiveOuterSteps`,
  default 5) keeps M1 dormant through the transient at tol 0.1, and the
  divergent setting then completes cleanly. The divergence was not re-run with
  the released build, but the hysteresis concerns when the reduction fires, so
  it remains required.
- **The flow gate limits firing on sloshing free surfaces.** A persistently
  sloshing free surface keeps `flowChange ~0.1-0.3` while the solid thermal
  field settles to ~1e-3, so the flow gate, not the thermal gate, decides. M1
  is therefore most effective for settled-interface or single-phase CHT; for
  sloshing VOF it trades a small accuracy cost for the speed-up under a loose
  flow tolerance (0.5 here, which also fires during parts of the start-up
  transient).
- **The first `flowChange` is huge** (e.g. 1.4e+15) because the initial
  velocity field is zero. This is harmless; it only keeps the first steps fully
  coupled.

### Usage recommendations

- Leave M1 off (the default) unless needed. For quasi-steady or long
  thermal-design runs, enable it with a conservative tolerance and keep
  `adaptiveOuterSteps >= 5`.
- A loose tolerance (0.5) gave -46% on the baseline settings and a further -40%
  on the preset at a few-hundredths of a kelvin QoI drift in this window; check
  the drift against an M1-off run before relying on it.
- When combining M1 with the preset (`nCorrectors 1`), monitor `alpha.water`.

## Regression

With M1 off the solver reproduces the previous release exactly: a 35-step
(0.005 s) run gives an identical solver log line by line (after filtering only
timing, run-header, dictionary-echo and `adaptiveOuter` lines), and the full
0.02 s window reproduces every iteration count, QoI and alpha bound of the
baseline (2776 `p_rgh` iterations, chip Tmax 302.9849 K), confirming that QW4
and the `fluid/pEqn.H` change are behaviour-neutral. The tuned preset on the
current solver reproduces the study's tuned result exactly in step count,
every iteration count, QoI and alpha bounds; only wall-clock differs.

## Production-scale levers (M5, M6)

Validated by partition/solver analysis, not by the serial window:

- **M6 (in preset).** The `fluid` region is ~95% of the cells and dominates the
  pressure-solve halo communication. Axis-aligned `hierarchical (4 2 2)` cut
  inter-processor faces 3.5% vs `scotch` at 16 ranks (8933 vs 9261, from
  `decomposePar -region fluid` on the 48x36x48 template); the gain grows with
  rank count. Re-tune `n` per rank count by comparing candidates (e.g.
  `(4 2 2)`, `(2 2 4)`, `(4 1 4)`) on the `Number of processor faces` reported
  by `decomposePar -region fluid`. The shipped `system/fluid/decomposeParDict`
  keeps `scotch`. Keep the **solid** regions on `scotch` at the same
  `numberOfSubdomains` - never under-decompose them, or the run FATALs at
  start-up.
- **M5 (production-gated, not in preset).** The serial test changed only
  `nCellsInCoarsestLevel` 10 -> 200 on `p_rgh` (2357 vs 2336 iters, +2%
  wall-clock; its `nFinestSweeps 2; nPostSweeps 2;` lines are the v2412 GAMG
  defaults): neutral-to-worse at 83k serial, with no parallel all-reduce to
  attack. A candidate production block to benchmark on the 663k / 16-rank
  mesh, on wall-clock rather than iterations (a recommendation only, not run in
  this study):

  ```
  p_rgh
  {
      ...
      nCellsInCoarsestLevel  500;
      processorAgglomerator  masterCoarsest;
  }
  ```

  In v2412 `processorAgglomerator` is only constructed when running on more
  than one processor, so it has no effect in serial.

## LTS and solid thermal preconditioning

- **LTS (declined).** Local time stepping accelerates a run to a *steady
  state*, but the M1 instrumentation shows `flowChange` never settles
  (~0.1-0.3 indefinitely): the impinging jet keeps the free surface oscillating,
  so this case has **no true steady state**. LTS would only suit a
  settled-interface or single-phase CHT configuration; the solver remains fully
  transient (`LTS` is a compile-time `false` in the driver).
- **Solid thermal preconditioning (recommended for steady-design runs).**
  Warm-starting the solids near their operating point (demonstrated stable with
  chip 322 K / heat-sink 312 K initial fields, alpha bounded) skips the slow
  conduction warm-up. The converged state is IC-independent, so this cuts
  steps-to-steady at zero numerical risk when the goal is the steady chip
  temperature rather than the cold-start transient.

## Tuned preset (opt-in)

The shipped example has the M1 keys present but disabled and the QW3
convergence line present but commented out in the top-level
`system/fvSolution`. The tuned preset validated in this study is the following
set of edits to the example; apply them yourself when you want the speed-up
(example 19 of
[openfoam-tutorial-ko](https://github.com/swtbkim/openfoam-tutorial-ko)
explains each one, in Korean).

| item | file | baseline (shipped) | tuned preset |
|---|---|---|---|
| QW2 | `system/fluid/fvSolution`, `solvers/p_rgh` | `relTol 0.01;` | `relTol 0.05;` (`p_rghFinal` stays `relTol 0`) |
| free | `system/fluid/fvSolution`, `"alpha.water.*"` | `nAlphaSubCycles 2;` | `nAlphaSubCycles 1;` |
| QW5a | `system/fluid/fvSchemes` | `div(rhoPhi,T) Gauss upwind;` | `div(rhoPhi,T) Gauss limitedLinear 1;` |
| QW3 | `system/fvSolution`, `PIMPLE/energyCoupling` | `iterations 2;` | `iterations 2; convergence { "T" 1e-3; }` |
| M2 | `system/fluid/fvSolution`, `PIMPLE` | `nCorrectors 2;` | `nCorrectors 1;` (see caveat) |
| M6 | `system/fluid/decomposeParDict` | `method scotch;` | `method hierarchical; coeffs { n (4 2 2); }` (16 ranks) |
| M1 | `system/fvSolution`, `PIMPLE` | `adaptiveOuterTol 0;` | `adaptiveOuterTol 0;` (off in the preset; optional `0.5` on top of it, see [M1](#m1---adaptive-outer-corrector-reduction-off-by-default)) |

The tuned run also set `cacheAgglomeration true;` on `p_rgh` explicitly; this
is the GAMG default and has no effect. Not part of the preset: QW1
(DICGaussSeidel) and M5, production-scale candidates to benchmark on the target
parallel run.

Net: **-27.8% wall-clock, accuracy-neutral, interface-bounded** with the preset
(-28.3% in the same-session M1 comparison); M1 on top of it gives -57% vs
baseline (see [Results](#results)). The production GAMG/decomposition levers
add headroom at 663k cells / 16 ranks.
