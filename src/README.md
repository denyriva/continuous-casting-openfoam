# Continuous Casting in OpenFOAM

Custom OpenFOAM Foundation v14 extensions for continuous casting with solidification, strand withdrawal, macrosegregation, composition-dependent phase equilibrium, and prescribed electromagnetic stirring (EMS).

This repository contains the custom source code developed for a staged validation project whose objective is to build an OpenFOAM continuous-casting model with capabilities approaching those normally available in commercial CFD packages, while keeping each added physical mechanism independently testable.

> **Current status:** Gate 5A completed and frozen as tag `gate5a-complete`.  
> Development of solutal buoyancy is being carried out on branch `gate5b-solutal-buoyancy`.

---

## 1. Scope

The model currently includes:

- incompressible melt flow;
- thermal buoyancy using a Boussinesq approximation;
- SST k-omega turbulence where required by the case;
- enthalpy-porosity solidification;
- mushy-zone momentum damping;
- latent heat;
- prescribed strand withdrawal / casting velocity;
- conservative treatment of latent-energy transport in the stationary casting domain;
- transport of a mixture carbon field;
- liquid/solid carbon partitioning;
- reconstruction of `CarbonL` and `CarbonS`;
- composition-dependent liquidus temperature;
- composition-dependent liquid fraction;
- macrosegregation caused by advection and solidification;
- optional prescribed Lorentz-force EMS source.

The present alloy implementation is a **binary Fe-C-style model** using carbon as the transported solute.

---

## 2. Repository structure

```text
continuous-casting-openfoam/
├── README.md
├── .gitignore
└── src/
    ├── castingSolidificationMelting/
    │   ├── castingSolidificationMelting.C
    │   ├── castingSolidificationMelting.H
    │   └── Make/
    │       ├── files
    │       └── options
    │
    └── macrosegregationTransport/
        ├── macrosegregationTransport.C
        ├── macrosegregationTransport.H
        └── Make/
            ├── files
            └── options
```

Compiled objects and platform-specific `wmake` output are intentionally excluded from version control.

---

## 3. OpenFOAM version

Development and validation have been performed using:

- **OpenFOAM Foundation v14**
- Ubuntu 24.04 under WSL
- GCC toolchain supplied with the OpenFOAM environment

The code should therefore be considered written specifically against the **Foundation v14 API** unless adapted for another release.

Before compiling, load the OpenFOAM environment in the usual way and verify:

```bash
foamVersion
```

---

# 4. Custom models

## 4.1 `castingSolidificationMelting`

`castingSolidificationMelting` is a custom `fvModel` derived from the concepts used by the Foundation `solidificationMelting` model and extended for continuous casting.

Its main responsibilities are:

1. determining liquid fraction;
2. applying mushy-zone resistance;
3. accounting for latent heat;
4. applying thermal buoyancy;
5. representing strand withdrawal;
6. supporting composition-dependent liquidus behavior;
7. exposing fields required by the macrosegregation model.

Important generated or maintained fields include:

```text
solidificationMelting1:alpha1
Teq
CarbonL
CarbonS
```

where:

- `alpha1` is the liquid fraction;
- `Teq` is the local composition-dependent equilibrium temperature used by the model;
- `CarbonL` is liquid-phase carbon concentration;
- `CarbonS` is solid-phase carbon concentration.

### Casting velocity

The strand withdrawal velocity is prescribed through the project casting-kinematics configuration and used consistently by the solidification and macrosegregation models.

The model does **not** physically translate the mesh. Instead, casting motion is represented in the governing equations.

### Conservative latent-energy treatment

The casting version uses a conservative latent contribution of the form

```text
-L [ d(rho*alpha)/dt + div(phi*alpha_f) ]
```

so that latent energy is treated consistently when solidifying material is advected through a stationary computational domain.

This is one of the main differences between this model and a simple stationary melting/solidification formulation.

---

## 4.2 `macrosegregationTransport`

`macrosegregationTransport` solves the transported mixture carbon field:

```text
Carbon
```

The project deliberately does **not** use the field name `C`, because OpenFOAM commonly uses `C` for cell-centre information.

The phase concentrations are reconstructed as:

```text
CarbonL
CarbonS
```

using the liquid fraction and the partition coefficient.

For a local liquid fraction `fl` and partition coefficient `kp`,

```text
CarbonS = kp * CarbonL
```

and the conserved mixture concentration satisfies

```text
Carbon = fl*CarbonL + (1 - fl)*CarbonS
```

therefore

```text
CarbonL = Carbon / [fl + kp*(1 - fl)]
```

The present formulation is therefore closest to a **local-equilibrium / Lever-rule limit** for phase partitioning.

---

# 5. Composition-dependent phase equilibrium

Gate 5A introduced two-way coupling between composition and solidification.

Previously the model effectively followed:

```text
temperature
    ↓
liquid fraction
```

while carbon was transported alongside the thermal solution.

The present model contains the additional path:

```text
Carbon
  ↓
CarbonL
  ↓
composition-dependent liquidus temperature
  ↓
liquid fraction
```

This closes an important macrosegregation feedback loop:

```text
Carbon
  ↓
CarbonL
  ↓
Tliq(CarbonL)
  ↓
liquid fraction
  ↓
mushy resistance / flow
  ↓
Carbon transport
```

The liquidus relationship currently uses a configurable linearized slope:

```text
Tliq = TliqRef + liquidusSlope*(CarbonL - CarbonRef)
```

For validation, deliberately exaggerated slopes may be used to make the coupling easy to verify. Production values should be based on appropriate alloy thermodynamic data.

---

# 6. Alloy-property configuration

The validation development uses a common alloy-property dictionary to avoid duplicating composition-dependent constants across different OpenFOAM dictionaries.

A representative configuration is:

```text
CarbonRef              0.01;
partitionCoefficient   0.25;

TsolRef                1593.15;
TliqRef                1688.15;

liquidusSlope          -1000;
```

The numerical values above are examples from validation work and **must not automatically be interpreted as production thermodynamic data**.

In particular, validation cases may intentionally use artificial parameter values so that coupling mechanisms produce clearly observable responses.

---

# 7. Building the libraries

The Git repository is now the authoritative source tree.

Do not edit older copies of the custom libraries elsewhere in the OpenFOAM directory tree.

## `castingSolidificationMelting`

```bash
cd ~/OpenFOAM/denyr-14/projects/openfoam-continuous-casting/src/castingSolidificationMelting
wmake libso
```

## `macrosegregationTransport`

```bash
cd ~/OpenFOAM/denyr-14/projects/openfoam-continuous-casting/src/macrosegregationTransport
wmake libso
```

The resulting shared libraries are written to the normal OpenFOAM user-library location, typically:

```bash
$FOAM_USER_LIBBIN
```

Check with:

```bash
echo $FOAM_USER_LIBBIN
```

and verify the libraries exist:

```bash
ls $FOAM_USER_LIBBIN | grep -E 'castingSolidificationMelting|macrosegregationTransport'
```

Expected library names are of the form:

```text
libcastingSolidificationMelting.so
libmacrosegregationTransport.so
```

---

# 8. Loading the libraries in a case

A case using the custom models must load the required libraries in `system/controlDict`, for example:

```text
libs
(
    "libcastingSolidificationMelting.so"
    "libmacrosegregationTransport.so"
);
```

Exact usage depends on the validation or production case.

---

# 9. Carbon field convention

The conserved composition field is:

```text
Carbon
```

A typical initial field is therefore placed in:

```text
0/Carbon
```

with an internal field such as:

```text
internalField uniform 0.01;
```

Boundary conditions must be chosen consistently with the physical case.

For casting cases, the nominal inlet composition should normally be consistent with the configured `CarbonRef`.

---

# 10. Numerical considerations

## Carbon transport

The carbon equation requires appropriate discretization and solver entries.

A representative convection scheme is:

```text
div(phiCarbon,Carbon) Gauss upwind;
```

and the scalar requires a solver entry in `fvSolution`.

The exact schemes used in final simulations should be chosen according to mesh quality, stability and accuracy requirements.

## Time step

Casting, buoyancy, EMS and segregation can produce strongly coupled transients.

Use a time step compatible with the velocity Courant number and, where appropriate, `adjustableRunTime` / `maxCo`.

## Mushy-zone damping

The solidification model uses a Darcy-like momentum resistance in partially solid cells.

The corresponding coefficient strongly affects how rapidly velocity is suppressed as liquid fraction falls and should therefore be treated as a physical/numerical model parameter rather than an arbitrary tuning constant.

---

# 11. Validation history

Development has been intentionally divided into gates so that each added mechanism is validated before being introduced into the full continuous-casting case.

## Gate 1 — Solidification sanity checks

Validated basic phase-change behavior and mushy-zone evolution.

## Gate 2 — Continuous casting

Added:

- mapped thermal initialization;
- strand withdrawal;
- conservative latent-energy transport;
- thermal buoyancy;
- three-dimensional casting flow;
- turbulence.

The developed no-EMS casting solution became the baseline for later work.

## Gate 3 — Electromagnetic stirring

Validated the prescribed EMS source independently before coupling it to the caster.

Checks included:

- force magnitude;
- direction/sign;
- spatial distribution;
- net force;
- torque;
- resulting flow rotation;
- interaction with thermal buoyancy.

The present EMS implementation is a **prescribed Lorentz-force model**. It does not solve Maxwell's equations internally.

## Gate 4 — Macrosegregation

Added:

- `Carbon`;
- `CarbonL`;
- `CarbonS`;
- partitioning;
- carbon transport;
- macrosegregation in the full caster.

No-EMS and EMS cases were compared after temporal development.

The tested EMS intensity did not improve segregation in the studied configuration and produced stronger segregation in the comparison performed during validation.

## Gate 5A — Composition-dependent phase equilibrium

Added:

- composition-dependent liquidus behavior;
- coupling of `CarbonL` to phase equilibrium;
- corresponding composition-dependent liquid fraction;
- dedicated 1D / no-flow validation.

Gate 5A is frozen in Git as:

```text
gate5a-complete
```

---

# 12. Current development: Gate 5B

Gate 5B is intended to add **solutal buoyancy**.

The planned extension is to modify the Boussinesq density contribution from a purely thermal form to a thermosolutal form, conceptually:

```text
rho ≈ rhoRef * [1 - betaT*(T - TRef) - betaC*(CarbonL - CarbonRef)]
```

This introduces the additional feedback:

```text
CarbonL
  ↓
liquid density
  ↓
solutal buoyancy
  ↓
flow
  ↓
Carbon transport
```

Development is isolated on:

```text
gate5b-solutal-buoyancy
```

The intended validation sequence is:

1. solutal buoyancy in a simple liquid domain;
2. combined thermal + solutal buoyancy;
3. full continuous-casting coupling.

---

# 13. Physics not yet included

The current model should not be interpreted as a complete industrial steel-solidification package.

Important effects not presently included include:

- finite solid-state back diffusion;
- Scheil-type microsegregation;
- multicomponent alloy chemistry;
- variable partition coefficients from full thermodynamic databases;
- full CALPHAD coupling;
- grain nucleation and grain-structure evolution;
- equiaxed crystal transport;
- shrinkage-driven flow;
- porosity formation;
- free-surface / mold-powder dynamics;
- air-gap formation;
- contact-resistance evolution;
- deformation of the solid shell;
- internally solved electromagnetic fields;
- induced-current calculation;
- Joule heating.

These are possible future extensions, not requirements for using the current binary macrosegregation model.

---

# 14. EMS model limitation

EMS is presently represented as a prescribed Lorentz-force field:

```text
FL = FL(x,t)
```

which is inserted into the momentum equation.

This is suitable when electromagnetic forces are known analytically or obtained from an external electromagnetic model.

The implementation does **not** currently solve the electromagnetic field equations or their mutual coupling with the fluid.

Therefore the EMS implementation should be understood as a reduced-order MHD coupling rather than a complete electromagnetic solver.

---

# 15. Git workflow

The GitHub repository is the source of truth for the custom code.

Current important references:

```text
main
gate5a-complete
gate5b-solutal-buoyancy
```

Recommended development workflow:

```bash
git status
git diff

# edit / compile / validate

git add .
git commit -m "Describe the implemented change"
git push
```

Before starting a new major model extension, create a dedicated branch or tag a known working state.

Example:

```bash
git tag some-stable-milestone
git push origin some-stable-milestone
```

This keeps validated physics reproducible while allowing experimental development to proceed safely.

---

# 16. Reproducing the development environment

A minimal setup workflow is:

```bash
# Clone
git clone git@github.com:denyriva/continuous-casting-openfoam.git

cd continuous-casting-openfoam

# Build solidification model
cd src/castingSolidificationMelting
wmake libso

# Build macrosegregation transport
cd ../macrosegregationTransport
wmake libso
```

A simulation case must then provide the appropriate:

```text
0/
constant/
system/
```

dictionaries and load the custom libraries.

Validation cases are not yet part of this repository.

---

# 17. Project philosophy

The project follows a deliberately staged development approach:

> **Do not add a physical coupling to the full caster until that coupling works in isolation.**

Each major mechanism is therefore tested in progressively more complex configurations before being accepted into the complete continuous-casting model.

This has been particularly important for:

- strand withdrawal;
- buoyancy;
- EMS;
- carbon transport;
- phase partitioning;
- composition-dependent liquidus behavior.

The same strategy is intended for future solutal-buoyancy and microsegregation extensions.

---

# 18. Disclaimer

This repository is a research/development implementation.

Validation cases have been designed primarily to verify mathematical implementation, coupling direction, conservation and qualitative physical behavior. Some test parameters are intentionally artificial.

Before using the model for engineering prediction, alloy-specific thermophysical and thermodynamic properties, numerical sensitivity, mesh convergence and comparison against experimental or industrial data should be established.

---

# 19. License and upstream code

This project contains custom code developed against OpenFOAM Foundation v14 and may include code structurally derived from or adapted from OpenFOAM components.

Before publishing the repository publicly, the applicable OpenFOAM licensing requirements and attribution should be reviewed and an appropriate repository license should be added.

OpenFOAM Foundation source code is distributed under the GNU General Public License.

---

## Status summary

| Capability | Status |
|---|---|
| Thermal solidification | Implemented |
| Mushy-zone damping | Implemented |
| Strand withdrawal | Implemented |
| Conservative latent-energy transport | Implemented |
| Thermal buoyancy | Implemented |
| SST turbulence | Supported |
| Carbon transport | Implemented |
| Liquid/solid partitioning | Implemented |
| Macrosegregation | Implemented |
| Composition-dependent liquidus | Implemented |
| Composition-dependent liquid fraction | Implemented |
| Prescribed EMS | Implemented |
| Solutal buoyancy | In development |
| Finite back diffusion | Not implemented |
| Multicomponent alloy | Not implemented |
| Full electromagnetic solution | Not implemented |
