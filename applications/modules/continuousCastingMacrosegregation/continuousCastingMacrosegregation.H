/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derived from the OpenFOAM Foundation v14 incompressibleFluid
    solver module and is distributed under the GNU General Public License.

Class
    Foam::solvers::continuousCastingMacrosegregation

Description
    Incompressible binary-alloy macrosegregation solver module.

    Current implementation stage:
      - standard incompressible pressure/velocity infrastructure;
      - alloy properties read from constant/alloyProperties;
      - prescribed solidVelocity read from constant/alloyProperties;
      - CC-2 moving-solid BKC momentum resistance based on U-solidVelocity;
      - temperature transport with dictionary-driven phase properties and
        phase-fraction-weighted Cp and k, with Eq. (13) sensible-energy
        advection discretized on the transported product Cp*T;
      - v11 adds a diagnostic-only audit of the exact discrete energy operator
        used by the final solidification correction;
      - Lever-rule phase equilibrium fields fs, CarbonL and CarbonS;
      - analytical Lever-rule dfs/dT retained as a diagnostic field;
      - direct discrete latent storage L*(fs-fs.oldTime())/dt with
        configurable fixed-point solidification iterations;
      - thermo-solutal Boussinesq buoyancy in momentum;
      - standard Blake-Kozeny-Carman (BKC) mushy-zone permeability and
        implicit Darcy damping toward the prescribed solid velocity;
      - Bennon-Incropera mixture species transport, discretized term-by-term
        as written in Eq. (19), with moving-solid relative advection based on
        U-us;
      - CC-8 turbulence-consistent heat and species transport using the
        momentum-model turbulent kinematic viscosity, with Prt and Sct read
        from constant/alloyProperties. Turbulent scalar transport is applied
        only to the liquid fraction.

    Hybrid-model terms are added in later validation stages.

\*---------------------------------------------------------------------------*/

#ifndef continuousCastingMacrosegregation_H
#define continuousCastingMacrosegregation_H

#include "basicFluidSolver.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "pressureReference.H"
#include "MRFZones.H"
#include "IOdictionary.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                 Class continuousCastingMacrosegregation Declaration
\*---------------------------------------------------------------------------*/

class continuousCastingMacrosegregation
:
    public basicFluidSolver
{

protected:

    // Pressure

        //- Pressure field
        volScalarField p_;

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Kinematic properties

        //- Velocity field
        volVectorField U_;

        //- Volumetric-flux field
        surfaceScalarField phi_;


    // Alloy properties

        //- Binary-alloy material-property dictionary
        IOdictionary alloyProperties_;

        //- Optional runtime controls for custom diagnostic output
        //  read from constant/diagnosticsProperties
        IOdictionary diagnosticsProperties_;

        //- Experimental pseudo-steady continuation controls
        IOdictionary pseudoSteadyProperties_;
        bool pseudoSteadyEnabled_;
        scalar pseudoMinDeltaT_;
        scalar pseudoMaxDeltaT_;
        scalar pseudoMaxCo_;
        scalar pseudoGrowthFactor_;
        scalar pseudoShrinkFactor_;
        scalar pseudoGrowthThreshold_;
        scalar pseudoMaxDeltaTStep_;
        scalar pseudoMaxDeltaCarbonStep_;
        scalar pseudoMaxDeltaFsStep_;
        scalar pseudoMaxDeltaUStep_;
        mutable scalar pseudoDeltaTTarget_;
        label pseudoReportInterval_;

        //- Reference/constant density [kg/m3]
        scalar rho_;

        //- Liquid specific heat capacity [J/(kg K)]
        scalar CpLiquid_;

        //- Solid specific heat capacity [J/(kg K)]
        scalar CpSolid_;

        //- Liquid thermal conductivity [W/(m K)]
        scalar kLiquid_;

        //- Solid thermal conductivity [W/(m K)]
        scalar kSolid_;

        //- Latent heat of fusion [J/kg]
        scalar latentHeat_;

        //- Initial/reference mixture solute mass fraction
        scalar Carbon0_;

        //- Equilibrium partition coefficient Cs/Cl
        scalar kp_;

        //- Liquidus slope [K per unit solute mass fraction]
        scalar liquidusSlope_;

        //- Melting temperature of the pure solvent [K]
        scalar Tmelt_;

        //- Liquid dynamic viscosity [Pa s]
        scalar muLiquid_;

        //- Liquid solute diffusivity [m2/s]
        scalar DL_;

        //- Solid solute diffusivity [m2/s]
        scalar DS_;

        //- Turbulent Prandtl number for liquid-phase heat transport
        scalar Prt_;

        //- Turbulent Schmidt number for liquid-phase species transport
        scalar Sct_;

        //- Thermal expansion coefficient [1/K]
        scalar betaT_;

        //- Solutal expansion coefficient [per unit solute mass fraction]
        scalar betaC_;

        //- Reference temperature used by the Boussinesq model [K]
        scalar TRef_;

        //- Secondary dendrite-arm spacing [m]
        scalar lambda2_;

        //- Microsegregation closure: lever or vollerBeckermann
        word microsegregationModel_;

        //- Voller-Beckermann coarsening correction alpha_C [-]
        scalar vollerBeckermannAlphaC_;

        //- Small regularisation for the vanishing-solid limit [-]
        scalar microsegregationXi_;

        //- Prescribed solid velocity [m/s]
        //  CC-2: used by moving-solid BKC momentum resistance
        vector solidVelocity_;

        //- Number of fixed-point solidification/latent-heat coupling loops
        label nSolidificationLoops_;

        //- Gravitational acceleration read from constant/g
        uniformDimensionedVectorField g_;


    // Thermal field

        //- Temperature field [K]
        volScalarField T_;


    // Binary-alloy phase/composition fields

        //- Mixture solute mass fraction
        volScalarField Carbon_;

        //- Solid mass fraction
        volScalarField fs_;

        //- Liquid-phase solute mass fraction
        volScalarField CarbonL_;

        //- Solid-phase average solute mass fraction
        volScalarField CarbonS_;

        //- Interfacial solid concentration C_s^* = kp*C_l^*
        volScalarField CarbonSInterface_;

        //- Local solidification time used by the Voller-Beckermann
        //  solutal Fourier number [s]
        volScalarField localSolidificationTime_;

        //- Voller-Beckermann back-diffusion coefficient beta [-]
        volScalarField backDiffusionCoefficient_;

        //- Effective mixture solute diffusivity used by the implicit
        //  mixture-C diffusion operator:
        //
        //    Dmix = fs*DS + fl*(DL + nut/Sct)
        //
        //  The accompanying explicit liquid/solid correction operators keep
        //  the final diffusion flux phase-separated.
        //  [m2/s]
        volScalarField speciesDiffusivity_;

        //- Algebraic ratio CarbonL/Carbon.
        //  Retained as an AUTO_WRITE diagnostic for the Lever partitioning:
        //    CarbonL = liquidAdvectionFactor*Carbon
        //  v9 no longer uses this ratio to collapse the two Eq. (19)
        //  advection terms into a single div(U CarbonL) operator.
        volScalarField liquidAdvectionFactor_;

        //- Macrosegregation degree
        //    100*(Carbon-Carbon0)/Carbon0
        //  [%]
        volScalarField macrosegregation_;

        //- Mixture specific heat capacity
        //    Cp = fs*CpSolid + (1-fs)*CpLiquid
        //  [J/(kg K)]
        volScalarField CpMix_;

        //- Effective thermal conductivity. For the AFRODITE benchmark
        //  rho_s = rho_l, hence gs = fs and
        //    kEff = fs*kSolid + (1-fs)*kLiquid
        //  [W/(m K)]
        volScalarField kEff_;

        //- Analytical Lever-rule derivative dfs/dT [1/K].
        //  Zero outside the mushy interval. Retained for diagnostics; the
        //  energy equation uses the direct discrete change in fs.
        volScalarField dfsdT_;

        //- Local apparent storage coefficient
        //    CpApp = CpMix - L*dfs/dT
        //  [J/(kg K)]
        //
        //  This is retained as a diagnostic field. The energy equation
        //  itself keeps the conservative d(CpMix*T)/dt term and uses the
        //  direct discrete latent source L*(fs-fs.oldTime())/dt.
        volScalarField CpApp_;


    // Thermo-solutal buoyancy fields

        //- Thermal Boussinesq acceleration [m/s2]
        volVectorField buoyancyThermal_;

        //- Solutal Boussinesq acceleration [m/s2]
        volVectorField buoyancySolutal_;

        //- Total Boussinesq acceleration [m/s2]
        volVectorField buoyancyTotal_;


    // Standard-BKC mushy-zone resistance fields

        //- Inverse permeability K^-1 [1/m2]
        volScalarField bkcInvK_;

        //- Kinematic Darcy damping coefficient nu_l/K [1/s]
        volScalarField bkcDragCoeff_;

        //- Darcy acceleration -nu_l/K * (U-solidVelocity) [m/s2]
        volVectorField bkcDragAcceleration_;


    // Momentum transport

        //- Kinematic viscosity model
        autoPtr<viscosityModel> viscosity;

        //- Pointer to the momentum transport model
        autoPtr<incompressible::momentumTransportModel> momentumTransport;


    // Optional models

        //- Reference to the MRFZones MeshObject
        const MRFZones& MRF;


    // Cached temporary fields

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> Uf;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;


    // Protected Member Functions

        //- Set rDeltaT for LTS
        void setRDeltaT();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Validate the alloy-property dictionary
        void validateAlloyProperties() const;

        //- Return true when the named diagnostic is enabled for the
        //  current time index according to constant/diagnosticsProperties
        bool diagnosticEnabled(const word& name) const;

        //- Validate experimental pseudo-steady controls
        void validatePseudoSteadyProperties() const;

        //- Update the target for the next pseudo time step
        void updatePseudoSteadyControl();

        //- Advance the local solidification-time history once per physical
        //  time step. Dong et al. neglect remelting, so the history is
        //  monotonic and freezes after complete solidification.
        void updateLocalSolidificationTime();

        //- Update phase state and microsegregation closure:
        //    Tliq = Tmelt + mliq*C
        //    Tsol = Tmelt + mliq*C/kp
        //    fs   = current composition-dependent phase relation
        //    lever: Cl=C/[1+fs(kp-1)], Cs=kp*Cl
        //    vollerBeckermann: Dong et al. Eqs. (17)-(23)
        //
        //  Returns the global maximum change in fs from the previous
        //  iteration.  If report=true, prints the full phase-state audit.
        scalar updatePhaseState(const bool report);

        //- Update thermo-solutal Boussinesq accelerations and diagnostics
        void updateBuoyancy();

        //- Update moving-solid BKC permeability / Darcy damping fields.
        //  If report=true, print the global BKC audit.
        void updateBKCDrag(const bool report);

        //- Update species-transport coefficient/diagnostic fields and,
        //  optionally, print the global solute-conservation audit.
        void updateSpeciesDiagnostics(const bool report);

        //- Solve the Bennon-Incropera mixture species equation with the
        //  moving-solid relative-advection term based on U-us.
        //  If report=true, print the CC-9 open-domain/discrete species audit
        //  immediately after the solve, before phase-state re-closure changes
        //  CarbonL/CarbonS or the frozen diffusion coefficients.
        //  Returns the global maximum change in Carbon during this correction.
        scalar solveSpeciesTransport(const bool report);

        //- Diagnose the relative contribution of temperature and mixture
        //  composition changes to Lever-rule solid-fraction evolution:
        //
        //    dfs/dt = (dfs/dT)_C dT/dt + (dfs/dC)_T dC/dt
        //
        //  Enabled explicitly with `phaseCoupling true;` in
        //  constant/diagnosticsProperties.
        void updatePhaseCouplingDiagnostics() const;

        //- Print a global Eq. (13) energy-conservation audit.
        //
        //  For standard BKC with us=0 and a closed cavity:
        //    d/dt integral[rho(Cp*T - L*fs)dV]
        //      = integral_boundary[k grad(T).n dA]
        void updateEnergyDiagnostics() const;

        //- Construct the pressure equation
        //  and correct the pressure and velocity
        virtual void correctPressure();


public:

    // Public Data

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the velocity field
        const volVectorField& U;

        //- Reference to the volumetric-flux field
        const surfaceScalarField& phi;

        //- Reference to the temperature field
        const volScalarField& T;

        //- Reference to the mixture solute field
        const volScalarField& Carbon;

        //- Reference to the solid mass-fraction field
        const volScalarField& fs;

        //- Reference to the liquid-phase solute field
        const volScalarField& CarbonL;

        //- Reference to the solid-phase average solute field
        const volScalarField& CarbonS;

        //- Reference to interfacial solid concentration
        const volScalarField& CarbonSInterface;

        //- Reference to local solidification time
        const volScalarField& localSolidificationTime;

        //- Reference to back-diffusion coefficient
        const volScalarField& backDiffusionCoefficient;

        //- Reference to the mixture solute diffusivity
        const volScalarField& speciesDiffusivity;

        //- Reference to CarbonL/Carbon
        const volScalarField& liquidAdvectionFactor;

        //- Reference to macrosegregation degree [%]
        const volScalarField& macrosegregation;

        //- Reference to the mixture specific-heat field
        const volScalarField& CpMix;

        //- Reference to the effective thermal-conductivity field
        const volScalarField& kEff;

        //- Reference to the analytical Lever-rule derivative dfs/dT
        const volScalarField& dfsdT;

        //- Reference to the apparent heat-capacity diagnostic field
        const volScalarField& CpApp;

        //- Reference to the inverse standard-BKC permeability
        const volScalarField& bkcInvK;

        //- Reference to the standard-BKC damping coefficient
        const volScalarField& bkcDragCoeff;

        //- Reference to the standard-BKC drag acceleration
        const volVectorField& bkcDragAcceleration;


    //- Runtime type information
    TypeName("continuousCastingMacrosegregation");


    // Constructors

        //- Construct from region mesh
        continuousCastingMacrosegregation(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        continuousCastingMacrosegregation
        (
            const continuousCastingMacrosegregation&
        ) = delete;


    //- Destructor
    virtual ~continuousCastingMacrosegregation();


    // Member Functions

        //- Return the current maximum time-step. In pseudo-steady mode the
        //  pseudoSteadyProperties Courant ceiling is used.
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Solve the temperature equation
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const continuousCastingMacrosegregation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
