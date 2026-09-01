# continuousCastingMacrosegregation

`continuousCastingMacrosegregation` is an OpenFOAM Foundation v14 solver module for transient macrosegregation in continuous casting of binary alloys. It couples incompressible flow, moving-solid Darcy resistance, thermo-solutal buoyancy, energy conservation with latent heat, mixture-species transport, and equilibrium/non-equilibrium microsegregation closures.

The current implementation is aimed at continuous-casting calculations in which the liquid, mushy and solid regions coexist and the casting strand moves with a prescribed solid velocity. It is derived from the Foundation `incompressibleFluid` solver-module architecture and is executed through `foamRun`.

## Current model capabilities

The present solver includes:

- transient incompressible momentum and pressure coupling through the Foundation PIMPLE framework;
- prescribed casting/solid velocity `solidVelocity`;
- moving-solid Blake-Kozeny-Carman (BKC) momentum resistance based on `U - solidVelocity`;
- thermal and solutal Boussinesq buoyancy;
- phase-dependent thermal properties;
- latent heat with local phase-change storage and moving-solid latent-energy advection;
- mixture carbon transport following the Bennon-Incropera / Dong formulation;
- liquid and solid molecular diffusion and turbulent liquid diffusivity;
- selectable species-diffusion linearisation;
- Lever-rule or Voller-Beckermann microsegregation;
- adaptive nonlinear relaxation for the coupled `T-C-fs` fixed-point iteration;
- adaptive **thermophysical physical-time subcycling** when the nonlinear coupling cannot converge over the full physical time step;
- extensive optional diagnostics for coupling, energy, species, buoyancy, drag, inventory and phase evolution.

The solver does **not** accept a physical time step whose thermophysical nonlinear coupling has failed. If the full thermophysical increment does not converge, only the thermophysical state is restored and the same physical interval is retried using smaller internal substeps.

---

## Numerical sequence of one physical time step

For every OpenFOAM physical time step `deltaT`, the main sequence is:

1. `preSolve()`
   - advance Voller-Beckermann physical solidification history from the previously accepted state;
   - update Courant number / local time-step support;
   - update the mesh if required.
2. Momentum transport prediction.
3. `momentumPredictor()`
   - update thermal and solutal buoyancy;
   - update moving-solid BKC resistance;
   - solve the momentum equation for `U`.
4. `thermophysicalPredictor()`
   - attempt the coupled temperature-carbon-solid-fraction problem over the full physical `deltaT`;
   - if necessary, restore the thermophysical checkpoint and retry using `2`, `4`, ... internal thermophysical substeps;
   - each successful internal substep must satisfy the nonlinear convergence criteria before the next internal substep is started.
5. `pressureCorrector()`
   - PIMPLE pressure corrections;
   - update `p`, `phi`, and `U`.
6. Momentum-transport correction
   - update viscosity / turbulence model.
7. `postSolve()`
   - refresh final BKC, flow, phase, species, macrosegregation and requested diagnostic quantities.

A detailed Mermaid flowchart is provided in [`logic.mmd`](logic.mmd).

---

## Momentum model

The momentum equation is written in kinematic-pressure form. Before assembly, the solver updates thermo-solutal buoyancy and the BKC drag coefficient.

### Moving-solid BKC resistance

The Darcy resistance acts on the velocity relative to the prescribed solid/casting velocity:

```text
a_D = -D (U - Us)
D   = nu_l / K
```

The permeability follows the standard BKC form used by the model:

```text
K^-1 = K0^-1 fs^2 / (1 - fs)^3
K0   = lambda2^2 / 180
```

where `lambda2` is the secondary dendrite-arm spacing.

The `-D U` contribution is implicit through `fvm::Sp`, while `+D Us` is an explicit source.

### Thermo-solutal buoyancy

The solver evaluates thermal and composition-dependent Boussinesq accelerations from the current thermophysical state. Thermal and solutal contributions are stored separately and summed in `buoyancyTotal` before the momentum solve.

---

## Energy equation

The current moving-solid energy formulation keeps the two latent-energy advection operators separate at the discrete level. After division by constant density, the modeled balance is conceptually

```text
d(CpMix T)/dt + div(U CpMix T)
    = div((kEff/rho) grad(T))
    + L d(fs)/dt
    + div(L fs U)
    - div(L fs (U - Us))
```

This separation is intentional. Although some terms can be combined algebraically in the continuous equation, bounded finite-volume interpolation can make the corresponding discrete operators non-equivalent.

Within each nonlinear correction:

- `CpMix`, `kEff`, and phase-dependent latent coefficients are lagged from the current phase state;
- `T` is solved implicitly;
- temperature nonlinear relaxation is applied;
- the phase/microsegregation closure is then recomputed.

For fully liquid cells, turbulent thermal diffusion uses `nut/Prt`. In the mush and solid regions the phase-dependent molecular conductivity is used as implemented in the solver.

---

## Species / macrosegregation equation

The mixture carbon equation follows the Bennon-Incropera / Dong form implemented as separate finite-volume operators:

```text
d(C)/dt + div(U C)
    = div(DmixEff grad(C))
    + div(fl DlEff grad(Cl - C))
    + div(fs Ds grad(Cs - C))
    - div((U - Us)(Cl - C))
```

with

```text
DlEff   = DL + nut/Sct
DmixEff = fs DS + fl DlEff
fl      = 1 - fs
```

The bulk mixture advection uses `U`; only the relative phase-advection term uses `U - Us`.

Two diffusion linearisations are available:

### `mixtureCorrection`

Historical formulation. The mixture diffusion is implicit and phase-difference corrections are explicit.

### `phaseLinearized`

The local phase closure is frozen during a Picard correction using

```text
Cl = ql C
Cs = qs C
```

so that the direct phase-separated diffusion flux is split into an implicit `grad(C)` contribution plus explicit gradients of the phase closure factors. This changes the numerical linearisation, not the underlying continuous species equation.

---

## Microsegregation models

Set in `constant/alloyProperties` with:

```foam
microsegregationModel lever;
```

or

```foam
microsegregationModel vollerBeckermann;
```

### Lever rule

The phase state is obtained algebraically from local temperature and mixture carbon concentration using the specified partition coefficient and liquidus relation.

### Voller-Beckermann

The model adds finite solid-state back-diffusion through the Voller-Beckermann closure. `localSolidificationTime` is a **physical-time history variable**. It is advanced once per accepted OpenFOAM physical step from the previously accepted phase state.

Thermophysical subcycling does not repeatedly advance `localSolidificationTime`; the internal retries are a numerical integration mechanism inside the already defined physical time interval.

---

## Nonlinear `T-C-fs` coupling

The thermophysical equations are solved by a Picard-type fixed-point iteration.

For each nonlinear correction:

1. save the current iteration fields;
2. solve the energy equation for `T`;
3. update phase equilibrium / microsegregation from the new `T`;
4. assemble and solve mixture carbon transport;
5. update phase equilibrium / microsegregation again from the new `C`;
6. evaluate the maximum changes in `T`, `Carbon`, and `fs`.

For `solidificationIterationMode converged`, the iteration stops only when all three criteria are satisfied after at least `minSolidificationIterations`:

```text
max |T(k+1) - T(k)|       <= temperatureCouplingTolerance
max |C(k+1) - C(k)|       <= carbonCouplingTolerance
max |fs(k+1) - fs(k)|     <= solidFractionCouplingTolerance
```

The current production settings are:

```foam
solidificationIterationMode          converged;
minSolidificationIterations          2;
maxSolidificationIterations          60;

temperatureCouplingTolerance         0.01;
carbonCouplingTolerance              1e-06;
solidFractionCouplingTolerance       1e-05;
```

`fixed` mode remains available for legacy calculations and simply executes `nSolidificationLoops` corrections without a convergence test.

---

## Adaptive nonlinear relaxation

The solver can adapt the three Picard relaxation factors when the nonlinear residual stalls.

A normalized worst-component residual is defined as

```text
R = max(
      dT / temperatureCouplingTolerance,
      dC / carbonCouplingTolerance,
      dfs / solidFractionCouplingTolerance
    )
```

If `R` fails to improve sufficiently for a prescribed number of consecutive iterations, all three nonlinear relaxation factors are reduced, down to `minimumNonlinearRelaxation`.

Recommended/current controls:

```foam
temperatureNonlinearRelaxation       0.5;
speciesNonlinearRelaxation           0.5;
solidFractionNonlinearRelaxation     0.5;

adaptiveNonlinearRelaxation          true;
minimumNonlinearRelaxation           0.125;
nonlinearRelaxationReductionFactor   0.5;
nonlinearResidualStallRatio          0.95;
nonlinearResidualBadIterations       2;
```

With these settings the typical sequence is `0.5 -> 0.25 -> 0.125` when the normalized residual stalls.

---

## Thermophysical physical-time subcycling

This is the current fail-safe for nonlinear thermophysical stiffness.

At the start of `thermophysicalPredictor()`, the solver checkpoints the thermophysical state. It first attempts the coupled problem using the full physical time step:

```text
Nsub = 1
DeltaTthermo = DeltaTphysical
```

If the nonlinear coupling reaches `maxSolidificationIterations` without satisfying all three convergence criteria, the attempted thermophysical state is discarded. The checkpoint is restored and the **same physical interval** is retried as multiple equal thermophysical substeps:

```text
Nsub = 2  ->  DeltaTthermo = DeltaTphysical / 2
Nsub = 4  ->  DeltaTthermo = DeltaTphysical / 4
Nsub = 8  ->  DeltaTthermo = DeltaTphysical / 8
```

Each internal substep must converge before the next one begins. The local substep histories of `T`, `Carbon`, `fs`, and `CpMix` are then advanced to the converged substep state.

The global OpenFOAM `Time` object is **not** advanced during these internal substeps. Momentum, pressure, fluxes and turbulence therefore remain at the same stage of the outer segregated physical-time solution.

Current controls:

```foam
thermophysicalSubcycling             true;
maxThermophysicalSubcycles           8;
thermophysicalSubcycleFactor         2;
```

If the highest allowed subcycle level still fails, the solver restores the thermophysical checkpoint and terminates with a `FatalError`. The unconverged physical state is therefore never accepted.

Typical log sequence:

```text
COUPLING_SUBCYCLE ... action=substep-failed ... nSub=1 ...
COUPLING_SUBCYCLE ... action=restore-and-retry failedNSub=1 retryNSub=2 ...
COUPLING_CONTROL  ... thermoSub=1/2 ... status=converged ...
COUPLING_CONTROL  ... thermoSub=2/2 ... status=converged ...
COUPLING_SUBCYCLE ... action=accepted nSub=2 ...
```

A useful filter is:

```bash
grep -E \
'COUPLING_SUBCYCLE.*(substep-failed|restore-and-retry|accepted nSub=[248])' \
log.foamRun
```

---

## Important `alloyProperties` entries

A representative SWRH82B configuration currently used with the model is:

```foam
rho                             7340;
CpLiquid                        650;
CpSolid                         650;
kLiquid                         33.5;
kSolid                          24.7;
latentHeat                      231637;

Carbon0                         0.0081;
kp                              0.34;
liquidusSlope                   -7800;
Tmelt                           1808;

muLiquid                        0.00461;
DL                              2.335361441e-07;
DS                              7.641829566e-10;
Prt                             0.9;
Sct                             1.0;

betaT                           2.0e-4;
betaC                           1.10;
TRef                            1758;

lambda2                         1e-4;
solidVelocity                   (0 0 0.0275);

microsegregationModel           vollerBeckermann;
speciesDiffusionForm            phaseLinearized;
vollerBeckermannAlphaC          0.1;
microsegregationXi              1e-12;

nSolidificationLoops            3;
solidificationIterationMode     converged;
minSolidificationIterations     2;
maxSolidificationIterations     60;

temperatureCouplingTolerance    0.01;
carbonCouplingTolerance         1e-06;
solidFractionCouplingTolerance  1e-05;

temperatureNonlinearRelaxation       0.5;
speciesNonlinearRelaxation           0.5;
solidFractionNonlinearRelaxation     0.5;
adaptiveNonlinearRelaxation          true;
minimumNonlinearRelaxation           0.125;
nonlinearRelaxationReductionFactor   0.5;
nonlinearResidualStallRatio          0.95;
nonlinearResidualBadIterations       2;

thermophysicalSubcycling             true;
maxThermophysicalSubcycles           8;
thermophysicalSubcycleFactor         2;
```

The values above document the current continuous-casting configuration; they are not universal material constants for every alloy or casting process.

---

## Main fields

Key solved or reconstructed fields include:

| Field | Meaning |
|---|---|
| `U` | mixture velocity |
| `p` | kinematic pressure |
| `phi` | mixture volumetric flux |
| `T` | temperature |
| `Carbon` | mixture carbon mass fraction |
| `fs` | solid fraction |
| `CarbonL` | liquid-phase carbon |
| `CarbonS` | solid-phase carbon |
| `CarbonSInterface` | interfacial solid composition used by the microsegregation closure |
| `CpMix` | phase-weighted heat capacity |
| `kEff` | phase-weighted thermal conductivity |
| `dfsdT` | phase-response diagnostic / latent linearisation field |
| `CpApp` | apparent heat-capacity diagnostic |
| `localSolidificationTime` | physical solidification-age history for Voller-Beckermann |
| `macrosegregation` | macrosegregation diagnostic |
| `bkcDragCoeff` | moving-solid Darcy damping coefficient |
| `buoyancyThermal` | thermal buoyancy acceleration |
| `buoyancySolutal` | solutal buoyancy acceleration |
| `buoyancyTotal` | total buoyancy acceleration |

Additional diagnostics are created and written according to the solver and `constant/diagnosticsProperties` configuration.

---

## Diagnostics

The solver contains optional diagnostics for the coupled thermophysical problem, energy balance, species conservation, phase evolution, buoyancy, drag and flow behavior.

Particularly useful nonlinear-coupling log tags are:

```text
COUPLING_DIAG
COUPLING_MAX_T
COUPLING_MAX_C
COUPLING_MAX_FS
COUPLING_ADAPT
COUPLING_CONTROL
COUPLING_SUBCYCLE
```

Use `constant/diagnosticsProperties` to control detailed diagnostic output where supported.

---

## Build

From the solver-module directory:

```bash
wclean
wmake
```

The library is installed in the user's OpenFOAM library directory and is loaded by `foamRun` through the case configuration.

Example parallel execution:

```bash
mpirun -np 12 foamRun -parallel > log.foamRun
```

---

## Source layout

```text
continuousCastingMacrosegregation/
├── continuousCastingMacrosegregation.C
├── continuousCastingMacrosegregation.H
├── momentumPredictor.C
├── correctPressure.C
├── moveMesh.C
├── setRDeltaT.C
├── logic.mmd
└── Make/
    ├── files
    └── options
```

`continuousCastingMacrosegregation.C/.H` contain the thermophysical, phase, species, diagnostics and solver-module infrastructure. Momentum and pressure assembly are split into the corresponding helper source files.

---

## Current numerical philosophy

The model uses a segregated physical-time solution with strongly coupled thermophysical fixed-point iterations inside each physical step. The objective is to keep the outer OpenFOAM time integration conventional while giving the stiff `T-C-fs` subsystem enough nonlinear and temporal resolution to converge reliably.

The hierarchy is therefore:

```text
physical OpenFOAM time step
    -> flow predictor
    -> thermophysical nonlinear coupling
         -> adaptive Picard relaxation
         -> if required: thermophysical physical-time subcycling
    -> pressure correction
    -> turbulence / transport correction
    -> diagnostics and write
```

The full physical `deltaT` remains the time increment seen by the outer flow solution. Internal thermophysical substeps always sum exactly to that same physical interval.
