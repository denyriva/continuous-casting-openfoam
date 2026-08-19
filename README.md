# Continuous Casting OpenFOAM — Binary Alloy Macrosegregation Solver

This repository contains a custom OpenFOAM Foundation v14 solver module for binary-alloy solidification and macrosegregation. The current main implementation is the `binaryAlloyMacrosegregation` solver, developed and validated against the standard BKC formulation used in the AFRODITE Sn–3 wt% Pb benchmark described by Moeinirad & Amani.

The older `castingSolidificationMelting`-based implementation is retained in the repository as a legacy development path and should be treated separately from the validated BKC solver.

---

## Current validated model

The validated solver is located at:

```text
applications/modules/binaryAlloyMacrosegregation/
```

Runtime solver name:

```text
binaryAlloyMacrosegregation
```

The validated development milestone is tagged:

```text
afrodite-bkc-validated
```

The solver reproduces the principal ingredients of the standard BKC benchmark formulation:

- incompressible laminar flow;
- thermal and solutal buoyancy;
- lever-rule phase equilibrium;
- Darcy resistance in the mushy zone using a secondary dendrite arm spacing;
- phase-weighted heat capacity and thermal conductivity;
- latent heat treatment consistent with the benchmark energy equation;
- binary-alloy species transport with phase partitioning;
- five nonlinear solidification corrections per time step;
- macrosegregation diagnostics;
- discrete energy-balance diagnostics.

The AFRODITE benchmark was used as the principal validation case before considering continuous-casting-specific extensions such as prescribed solid withdrawal velocity.

---

## Repository structure

```text
.
├── applications/
│   └── modules/
│       └── binaryAlloyMacrosegregation/
│           ├── Make/
│           │   ├── files
│           │   └── options
│           ├── binaryAlloyMacrosegregation.C
│           ├── binaryAlloyMacrosegregation.H
│           ├── correctPressure.C
│           ├── momentumPredictor.C
│           ├── moveMesh.C
│           ├── setRDeltaT.C
│           └── logic.mmd
├── src/
│   └── ...
└── README.md
```

The `src/` tree contains the earlier fvModel-based continuous-casting implementation. It is preserved for reference and later reconciliation with the benchmarked solver, but it is not the implementation used for the AFRODITE BKC validation.

---

## OpenFOAM version

The solver was developed for:

```text
OpenFOAM Foundation v14
```

The development environment used WSL Ubuntu with the Foundation distribution of OpenFOAM.

The module is intended for execution through `foamRun`.

---

## Building the solver

From the solver-module directory:

```bash
cd ~/OpenFOAM/$USER-14/projects/openfoam-continuous-casting/applications/modules/binaryAlloyMacrosegregation
wmake libso .
```

The resulting library is written to the user's OpenFOAM library directory, typically:

```text
$FOAM_USER_LIBBIN/libbinaryAlloyMacrosegregationSolver.so
```

Check that it exists with:

```bash
ls -l $FOAM_USER_LIBBIN/libbinaryAlloyMacrosegregationSolver.so
```

---

## Running a case

A case using the solver should load the custom library and select the runtime solver.

Typical execution:

```bash
foamRun
```

Parallel execution:

```bash
decomposePar
mpirun -np <N> foamRun -parallel
```

For example:

```bash
mpirun -np 12 foamRun -parallel > log.foamRun
```

The case must provide the usual incompressible-flow fields together with temperature and alloy composition.

Typical primary fields:

```text
0/U
0/p
0/T
0/Carbon
```

The solver additionally creates and writes derived alloy fields such as:

```text
fs
CarbonL
CarbonS
macrosegregation
```

---

# Physical model

## 1. Flow

The liquid/mushy flow is treated as incompressible.

The momentum equation includes:

- pressure;
- viscous diffusion;
- thermal buoyancy;
- solutal buoyancy;
- BKC Darcy resistance in the mushy zone.

The validated AFRODITE case used laminar Stokes momentum transport.

---

## 2. Thermo-solutal buoyancy

The benchmark Boussinesq source is implemented in the form

$$
\mathbf{a}_b
=
-\mathbf{g}\,\beta_T\,(T-T_{ref})
-\mathbf{g}\,\beta_C\,(C_l-C_0).
$$

where:

- $T$ is temperature;
- $T_{ref}$ is the reference temperature;
- $C_l$ is liquid composition;
- $C_0$ is nominal alloy composition;
- $\beta_T$ is the thermal expansion coefficient;
- $\beta_C$ is the solutal expansion coefficient.

The AFRODITE Sn–Pb benchmark uses a negative solutal expansion coefficient, making Pb-rich liquid denser and therefore gravitationally unstable in the expected direction.

---

## 3. Phase equilibrium

The binary phase diagram is represented using linear liquidus and solidus relations.

Liquidus:

$$
T_{liq}
=
T_{melt}
+
m_{liq} C
$$

Solidus:

$$
T_{sol}
=
T_{melt}
+
m_{liq}\frac{C}{k_p}
$$

where:

- $T_{melt}$ is the pure-solvent melting temperature;
- $m_{liq}$ is the liquidus slope;
- $C$ is local mixture composition;
- $k_p$ is the partition coefficient.

In the mushy interval, the solid fraction follows the lever-rule expression

$$
f_s
=
\frac{1}{1-k_p}
\frac{T-T_{liq}}{T-T_{melt}}.
$$

The phase compositions are reconstructed from the mixture composition:

$$
C_l
=
\frac{C}
{1+f_s(k_p-1)}
$$

and

$$
C_s
=
k_p C_l.
$$

The implementation enforces the mixture closure

$$
C=f_s C_s+(1-f_s)C_l.
$$

---

## 4. BKC mushy-zone resistance

The standard BKC permeability relation is used.

$$
K_0
=
\frac{\lambda_2^2}{180}
$$

with $\lambda_2$ the secondary dendrite arm spacing.

The inverse permeability is

$$
K^{-1}
=
K_0^{-1}
\frac{f_{s,B}^2}
{(1-f_{s,B})^3}
$$

with the bounded solid fraction

$$
f_{s,B}=\min(f_s,0.99).
$$

The momentum sink is

$$
S_u
=
-\frac{\mu_l}{K}
(\mathbf{u}-\mathbf{u}_s).
$$

For the validated standard BKC benchmark,

$$
\mathbf{u}_s=0.
$$

The resistance is added implicitly to the momentum equation.

This is important: the validated AFRODITE implementation is the stationary-solid BKC model. A nonzero solid withdrawal velocity for continuous casting is a later extension and is not part of the benchmark tag.

---

## 5. Thermal properties

Heat capacity is phase weighted:

$$
c_p
=
f_s c_{p,s}
+
(1-f_s)c_{p,l}.
$$

Thermal conductivity is likewise phase weighted:

$$
k_{eff}
=
f_s k_s
+
(1-f_s)k_l.
$$

For the equal-density benchmark formulation, the volume and mass phase fractions coincide.

---

## 6. Energy equation

The implementation follows the benchmark energy formulation.

For the standard BKC case with stationary solid,

$$
\mathbf{u}_s=0,
$$

the latent convective terms in the published formulation cancel exactly.

The local latent contribution follows the benchmark approximation

$$
\frac{\partial(\rho L f_s)}{\partial t}
=
\rho L
\frac{\partial f_s}{\partial T}
\frac{\partial T}{\partial t}.
$$

Accordingly, the implementation uses the analytical lever-rule derivative $\partial f_s/\partial T$.

The solver does **not** introduce an additional

$$
\frac{\partial f_s}{\partial C}
\frac{\partial C}{\partial t}
$$

latent source, because that term is not part of the benchmark formulation being reproduced.

The sensible energy transport uses the phase-weighted product $c_p T$. To preserve the intended product-advection form while retaining an implicit finite-volume solve, the implementation uses an implicit split plus deferred correction.

---

## 7. Species transport

The mixture composition field is named:

```text
Carbon
```

This naming convention is deliberate; `C` is avoided because OpenFOAM commonly uses `C` for cell-centre coordinates.

Associated derived phase fields are:

```text
CarbonL
CarbonS
```

The species equation is implemented term-by-term following the benchmark equal-density binary-alloy form rather than algebraically collapsing the advection to a single `div(U*Cl)` operator.

This distinction matters numerically: two continuously equivalent forms are not necessarily equivalent under the same finite-volume discretization.

The implemented structure contains:

- mixture composition transient term;
- mixture advection;
- mixture diffusivity;
- liquid-phase diffusion correction;
- solid-phase diffusion correction;
- relative liquid-composition advection correction.

The phase-weighted mixture diffusivity is

$$
D
=
f_s D_s
+
(1-f_s)D_l.
$$

For the AFRODITE benchmark:

$$
D_s=0.
$$

---

## 8. Macrosegregation

Macrosegregation is reported as

$$
MS
=
\frac{C-C_0}{C_0}\times 100\%.
$$

Positive values indicate local enrichment relative to the nominal alloy composition, and negative values indicate depletion.

---

# Nonlinear coupling

Solidification, thermal properties, latent heat, species partitioning, and buoyancy are strongly coupled.

The benchmark solver therefore performs multiple solidification corrections per time step.

The validated AFRODITE setup used:

```text
nSolidificationLoops 5;
```

Within these corrections the solver updates the phase state and associated thermophysical fields until the coupled temperature/composition/solid-fraction state is sufficiently converged for that time step.

---

# AFRODITE benchmark

The principal validation case is the AFRODITE Sn–3 wt% Pb cavity benchmark.

Geometry:

```text
100 mm × 60 mm × 10 mm
```

The numerical model is effectively two-dimensional, using one cell through thickness.

Validated mesh:

```text
100 × 60 × 1
```

for a total of:

```text
6000 cells
```

Typical gravity:

```text
(0 -9.81 0)
```

Initial alloy composition:

```text
C0 = 0.03
```

Initial temperature:

```text
531.75 K
```

The thermal boundary condition follows the AFRODITE benchmark heating/cooling schedule used in the reference implementation.

---

## AFRODITE material properties

The validated case used approximately:

```text
rho             7130 kg/m3

CpSolid          209 J/(kg K)
CpLiquid         243 J/(kg K)

kSolid            55 W/(m K)
kLiquid           33 W/(m K)

latentHeat     56140 J/kg

muLiquid        0.002 Pa s

DL              3.5e-9 m2/s
DS              0

betaT           9.5e-5 1/K
betaC          -0.53

lambda2         9.0e-5 m

Tmelt           505.15 K
liquidusSlope  -128.6089 K/(mass fraction)
kp              0.0656

TRef            531.75 K
Carbon0         0.03
```

These values belong to the AFRODITE benchmark configuration and should not be treated as generic steel properties.

---

# Recommended numerical settings

The validated benchmark configuration used:

```text
Time scheme:
    Euler

Gradient:
    Gauss linear

Interpolation:
    linear

Advection:
    bounded Gauss linearUpwind

Diffusion:
    Gauss linear corrected
```

Pressure:

```text
GAMG
tolerance 1e-6
```

Velocity, temperature and composition:

```text
PBiCGStab
DILU
tolerance 1e-8
```

Pressure-velocity coupling:

```text
PIMPLE/PISO-style
nCorrectors 2
```

Adaptive time stepping was used with approximately:

```text
maxCo       0.5
maxDeltaT   0.05
```

---

# Validation status

The current solver reproduces the principal behaviour of the standard BKC implementation reported for the AFRODITE benchmark.

At approximately 2250 s, the validated solution showed:

- a solidification-front position close to the published standard-BKC result;
- liquid velocities of the same order and profile as the benchmark;
- strong mushy-zone damping;
- the expected dense solute-rich lower region caused by Pb segregation;
- conserved global alloy inventory;
- converged nonlinear solidification corrections.

The remaining differences observed during validation were profile-level differences rather than evidence of a missing dominant model term.

---

## Energy-conservation validation

A dedicated discrete energy audit was added before freezing the validated BKC milestone.

The older field-based diagnostic showed an apparent global mismatch of roughly 5%, but this was traced to the diagnostic definition rather than the solved equation.

The discrete finite-volume Eq. (13) audit showed residuals on the order of approximately:

```text
~1e-3 W
```

against equation terms of approximately:

```text
20–25 W
```

corresponding to relative errors of only a few $10^{-5}$, i.e. of order:

```text
~0.004 %
```

on average in the final validation window.

This confirmed that the assembled energy equation is internally conservative to a level far below the benchmark-model uncertainty.

---

# Important implementation notes

## Discrete equivalence matters

During development, an algebraically simplified species-advection form was found to produce noticeably different hydrodynamics from the term-by-term benchmark equation when discretized with bounded `linearUpwind`.

The final solver therefore keeps the benchmark species equation in its explicit term-by-term form.

The same principle is used for the energy advection treatment.

---

## Stationary solid in the validated benchmark

The validated BKC milestone assumes:

```text
solid velocity = 0
```

This must not be confused with the intended continuous-casting extension, where the solid shell is withdrawn through the domain.

Any introduction of `pullVelocity` should be treated as a new model extension and validated independently against the frozen BKC baseline.

---

## Equal-density formulation

The AFRODITE benchmark used the equal-density binary-alloy simplification.

Therefore the current validated implementation should not automatically be interpreted as a full variable-density multicomponent solidification solver.

---

# Legacy implementation

The repository also contains the earlier custom fvModel development under:

```text
src/
```

That code includes the original continuous-casting-oriented `castingSolidificationMelting` model and associated macrosegregation developments.

It remains useful because it contains earlier work on:

- continuous-casting pull velocity;
- composition-dependent phase equilibrium;
- thermo-solutal buoyancy;
- electromagnetic stirring coupling;
- macrosegregation transport.

However, the new `binaryAlloyMacrosegregation` solver is the reference implementation for BKC physics and benchmark fidelity.

The legacy implementation should therefore be documented and maintained separately rather than mixed with the validated benchmark description.

---

# Development workflow

The validated BKC state is preserved by the Git tag:

```bash
git checkout afrodite-bkc-validated
```

Do not develop new continuous-casting features directly on the validation tag.

A recommended workflow for future development is:

```bash
git checkout main
git pull
git checkout -b <new-feature-branch>
```

Examples:

```text
continuous-casting-pull-velocity
moving-solid-bkc
continuous-casting-validation
```

This keeps the AFRODITE benchmark state recoverable at all times.

---

# Next development stage

The natural next step is to extend the validated BKC model toward continuous casting.

The principal addition will be a nonzero solid velocity:

$$
\mathbf{u}_s \neq 0.
$$

This affects more than the Darcy sink. It also changes the transport structure of the benchmark equations, particularly the relative momentum/species/latent transport terms.

For that reason, the continuous-casting extension should be implemented incrementally and checked against the frozen standard-BKC baseline after each change.

---

# Reference

The implementation was benchmarked against the standard BKC formulation described in:

**Moeinirad, S. & Amani, M. — “Systematic Benchmarking of Macrosegregation: The Performance of a Modified Hybrid Model.”**

The AFRODITE Sn–3 wt% Pb solidification benchmark described in that work was used as the principal reference case for the current solver.

---

# Status

Current project status:

```text
Standard BKC benchmark implementation: validated
AFRODITE benchmark: reproduced with good fidelity
Energy discrete balance: validated
Species inventory: conserved
Continuous-casting pull velocity: not yet added to the validated BKC solver
```

Validated Git milestone:

```text
afrodite-bkc-validated
```
