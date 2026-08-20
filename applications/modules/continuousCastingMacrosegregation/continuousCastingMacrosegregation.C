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

\*---------------------------------------------------------------------------*/

#include "continuousCastingMacrosegregation.H"

#include "localEulerDdtScheme.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(continuousCastingMacrosegregation, 0);
    addToRunTimeSelectionTable(solver, continuousCastingMacrosegregation, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::solvers::continuousCastingMacrosegregation::correctCoNum()
{
    basicFluidSolver::correctCoNum(phi);
}


void Foam::solvers::continuousCastingMacrosegregation::continuityErrors()
{
    basicFluidSolver::continuityErrors(phi);
}


void Foam::solvers::continuousCastingMacrosegregation::validateAlloyProperties() const
{
    if (rho_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "rho must be > 0. Current value: " << rho_
            << exit(FatalIOError);
    }

    if (CpLiquid_ <= SMALL || CpSolid_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "CpLiquid and CpSolid must both be > 0." << nl
            << "    CpLiquid = " << CpLiquid_ << nl
            << "    CpSolid  = " << CpSolid_
            << exit(FatalIOError);
    }

    if (kLiquid_ <= SMALL || kSolid_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "kLiquid and kSolid must both be > 0." << nl
            << "    kLiquid = " << kLiquid_ << nl
            << "    kSolid  = " << kSolid_
            << exit(FatalIOError);
    }

    if (latentHeat_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "latentHeat must be > 0. Current value: " << latentHeat_
            << exit(FatalIOError);
    }

    if (Carbon0_ < 0 || Carbon0_ > 1)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Carbon0 must satisfy 0 <= Carbon0 <= 1. Current value: "
            << Carbon0_ << exit(FatalIOError);
    }

    if (kp_ <= 0 || kp_ > 1)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "kp must satisfy 0 < kp <= 1. Current value: " << kp_
            << exit(FatalIOError);
    }

    if (mag(liquidusSlope_) <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "liquidusSlope must be non-zero."
            << exit(FatalIOError);
    }

    if (muLiquid_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "muLiquid must be > 0. Current value: " << muLiquid_
            << exit(FatalIOError);
    }

    if (DL_ < 0 || DS_ < 0)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "DL and DS must be non-negative. Current values: DL="
            << DL_ << ", DS=" << DS_
            << exit(FatalIOError);
    }

    if (lambda2_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "lambda2 must be > 0. Current value: " << lambda2_
            << exit(FatalIOError);
    }

    if (nSolidificationLoops_ < 1)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "nSolidificationLoops must be >= 1. Current value: "
            << nSolidificationLoops_
            << exit(FatalIOError);
    }
}


Foam::scalar
Foam::solvers::continuousCastingMacrosegregation::updatePhaseState
(
    const bool report
)
{
    scalar minT = GREAT;
    scalar maxT = -GREAT;

    scalar minCarbon = GREAT;
    scalar maxCarbon = -GREAT;

    scalar minTliq = GREAT;
    scalar maxTliq = -GREAT;

    scalar minTsol = GREAT;
    scalar maxTsol = -GREAT;

    scalar minFs = GREAT;
    scalar maxFs = -GREAT;

    scalar minCarbonL = GREAT;
    scalar maxCarbonL = -GREAT;

    scalar minCarbonS = GREAT;
    scalar maxCarbonS = -GREAT;

    scalar minCpMix = GREAT;
    scalar maxCpMix = -GREAT;

    scalar minKEff = GREAT;
    scalar maxKEff = -GREAT;

    scalar minDfsdT = GREAT;
    scalar maxDfsdT = -GREAT;

    scalar minCpApp = GREAT;
    scalar maxCpApp = -GREAT;

    scalar maxMassClosureError = 0.0;
    scalar maxDeltaFs = 0.0;

    label nLiquid = 0;
    label nMushy = 0;
    label nSolid = 0;

    forAll(T_, celli)
    {
        const scalar Tcell = T_[celli];
        const scalar Craw = Carbon_[celli];

        if (Craw < -SMALL || Craw > 1.0 + SMALL)
        {
            FatalErrorInFunction
                << "Invalid mixture solute mass fraction." << nl
                << "    cell   = " << celli << nl
                << "    Carbon = " << Craw << nl
                << "Expected 0 <= Carbon <= 1."
                << abort(FatalError);
        }

        // Remove only round-off excursions at the physical bounds.
        const scalar Ccell = min(max(Craw, scalar(0)), scalar(1));

        // Moeinirad & Amani, Eqs. (22) and (23)
        const scalar Tliq =
            Tmelt_ + liquidusSlope_*Ccell;

        const scalar Tsol =
            Tmelt_ + liquidusSlope_*Ccell/kp_;

        scalar fsCell = 0.0;
        scalar dfsdTCell = 0.0;

        if (Tcell >= Tliq)
        {
            fsCell = 0.0;
            dfsdTCell = 0.0;
        }
        else if (Tcell <= Tsol)
        {
            fsCell = 1.0;
            dfsdTCell = 0.0;
        }
        else
        {
            // Lever rule, Eq. (28)
            const scalar denominator =
                (1.0 - kp_)*(Tcell - Tmelt_);

            if (mag(denominator) <= SMALL)
            {
                FatalErrorInFunction
                    << "Lever-rule denominator is too small." << nl
                    << "    cell   = " << celli << nl
                    << "    T      = " << Tcell << " K" << nl
                    << "    Tmelt  = " << Tmelt_ << " K"
                    << abort(FatalError);
            }

            fsCell =
                (Tcell - Tliq)/denominator;

            fsCell =
                min(max(fsCell, scalar(0)), scalar(1));

            // Analytical derivative of Eq. (28) at fixed mixture
            // composition C.  This is the quantity used in the paper's
            // Eq. (16) for the latent-heat recast.
            //
            // dfs/dT = mliq*C / [(1-kp)*(T-Tmelt)^2]
            const scalar dTToMelt = Tcell - Tmelt_;

            dfsdTCell =
                liquidusSlope_*Ccell
               /((1.0 - kp_)*sqr(dTToMelt));
        }

        // Paper Eqs. (14) and (15). Since rho_s = rho_l for this
        // benchmark, solid mass fraction and solid volume fraction coincide.
        const scalar CpMixCell =
            fsCell*CpSolid_ + (1.0 - fsCell)*CpLiquid_;

        const scalar kEffCell =
            fsCell*kSolid_ + (1.0 - fsCell)*kLiquid_;

        // Diagnostic apparent storage coefficient associated with the
        // latent term recast of Eq. (16).
        const scalar CpAppCell =
            CpMixCell - latentHeat_*dfsdTCell;

        if (CpAppCell <= SMALL)
        {
            FatalErrorInFunction
                << "Non-positive apparent heat capacity." << nl
                << "    cell    = " << celli << nl
                << "    T       = " << Tcell << " K" << nl
                << "    Carbon  = " << Ccell << nl
                << "    dfsdT   = " << dfsdTCell << " 1/K" << nl
                << "    CpApp   = " << CpAppCell << " J/(kg K)"
                << abort(FatalError);
        }

        // Lever-rule phase compositions, Eqs. (29) and (30)
        const scalar compositionDenominator =
            1.0 + fsCell*(kp_ - 1.0);

        if (compositionDenominator <= SMALL)
        {
            FatalErrorInFunction
                << "Invalid Lever-rule composition denominator." << nl
                << "    cell        = " << celli << nl
                << "    fs          = " << fsCell << nl
                << "    kp          = " << kp_
                << abort(FatalError);
        }

        const scalar CarbonLCell =
            Ccell/compositionDenominator;

        const scalar CarbonSCell =
            kp_*CarbonLCell;

        const scalar fsOldIter = fs_[celli];

        fs_[celli] = fsCell;
        CarbonL_[celli] = CarbonLCell;
        CarbonS_[celli] = CarbonSCell;
        CpMix_[celli] = CpMixCell;
        kEff_[celli] = kEffCell;
        dfsdT_[celli] = dfsdTCell;
        CpApp_[celli] = CpAppCell;

        const scalar reconstructedCarbon =
            fsCell*CarbonSCell
          + (1.0 - fsCell)*CarbonLCell;

        const scalar massClosureError =
            mag(Ccell - reconstructedCarbon);

        maxDeltaFs =
            max(maxDeltaFs, mag(fsCell - fsOldIter));

        maxMassClosureError =
            max(maxMassClosureError, massClosureError);

        minT = min(minT, Tcell);
        maxT = max(maxT, Tcell);

        minCarbon = min(minCarbon, Ccell);
        maxCarbon = max(maxCarbon, Ccell);

        minTliq = min(minTliq, Tliq);
        maxTliq = max(maxTliq, Tliq);

        minTsol = min(minTsol, Tsol);
        maxTsol = max(maxTsol, Tsol);

        minFs = min(minFs, fsCell);
        maxFs = max(maxFs, fsCell);

        minCarbonL = min(minCarbonL, CarbonLCell);
        maxCarbonL = max(maxCarbonL, CarbonLCell);

        minCarbonS = min(minCarbonS, CarbonSCell);
        maxCarbonS = max(maxCarbonS, CarbonSCell);

        minCpMix = min(minCpMix, CpMixCell);
        maxCpMix = max(maxCpMix, CpMixCell);

        minKEff = min(minKEff, kEffCell);
        maxKEff = max(maxKEff, kEffCell);

        minDfsdT = min(minDfsdT, dfsdTCell);
        maxDfsdT = max(maxDfsdT, dfsdTCell);

        minCpApp = min(minCpApp, CpAppCell);
        maxCpApp = max(maxCpApp, CpAppCell);

        if (fsCell <= SMALL)
        {
            ++nLiquid;
        }
        else if (fsCell >= 1.0 - SMALL)
        {
            ++nSolid;
        }
        else
        {
            ++nMushy;
        }
    }

    fs_.correctBoundaryConditions();
    CarbonL_.correctBoundaryConditions();
    CarbonS_.correctBoundaryConditions();
    CpMix_.correctBoundaryConditions();
    kEff_.correctBoundaryConditions();
    dfsdT_.correctBoundaryConditions();
    CpApp_.correctBoundaryConditions();

    // Global reductions so the diagnostics remain valid in parallel.
    minT = returnReduce(minT, minOp<scalar>());
    maxT = returnReduce(maxT, maxOp<scalar>());

    minCarbon = returnReduce(minCarbon, minOp<scalar>());
    maxCarbon = returnReduce(maxCarbon, maxOp<scalar>());

    minTliq = returnReduce(minTliq, minOp<scalar>());
    maxTliq = returnReduce(maxTliq, maxOp<scalar>());

    minTsol = returnReduce(minTsol, minOp<scalar>());
    maxTsol = returnReduce(maxTsol, maxOp<scalar>());

    minFs = returnReduce(minFs, minOp<scalar>());
    maxFs = returnReduce(maxFs, maxOp<scalar>());

    minCarbonL = returnReduce(minCarbonL, minOp<scalar>());
    maxCarbonL = returnReduce(maxCarbonL, maxOp<scalar>());

    minCarbonS = returnReduce(minCarbonS, minOp<scalar>());
    maxCarbonS = returnReduce(maxCarbonS, maxOp<scalar>());

    minCpMix = returnReduce(minCpMix, minOp<scalar>());
    maxCpMix = returnReduce(maxCpMix, maxOp<scalar>());

    minKEff = returnReduce(minKEff, minOp<scalar>());
    maxKEff = returnReduce(maxKEff, maxOp<scalar>());

    minDfsdT = returnReduce(minDfsdT, minOp<scalar>());
    maxDfsdT = returnReduce(maxDfsdT, maxOp<scalar>());

    minCpApp = returnReduce(minCpApp, minOp<scalar>());
    maxCpApp = returnReduce(maxCpApp, maxOp<scalar>());

    maxMassClosureError =
        returnReduce(maxMassClosureError, maxOp<scalar>());

    maxDeltaFs =
        returnReduce(maxDeltaFs, maxOp<scalar>());

    nLiquid = returnReduce(nLiquid, sumOp<label>());
    nMushy = returnReduce(nMushy, sumOp<label>());
    nSolid = returnReduce(nSolid, sumOp<label>());

    if (report)
    {
        Info<< "Lever-rule phase-state audit" << nl
            << "    T min/max             = "
            << minT << " " << maxT << " K" << nl
            << "    Carbon min/max        = "
            << minCarbon << " " << maxCarbon << nl
            << "    Tliq min/max          = "
            << minTliq << " " << maxTliq << " K" << nl
            << "    Tsol min/max          = "
            << minTsol << " " << maxTsol << " K" << nl
            << "    fs min/max            = "
            << minFs << " " << maxFs << nl
            << "    CarbonL min/max       = "
            << minCarbonL << " " << maxCarbonL << nl
            << "    CarbonS min/max       = "
            << minCarbonS << " " << maxCarbonS << nl
            << "    CpMix min/max         = "
            << minCpMix << " " << maxCpMix << " J/(kg K)" << nl
            << "    kEff min/max          = "
            << minKEff << " " << maxKEff << " W/(m K)" << nl
            << "    dfsdT min/max         = "
            << minDfsdT << " " << maxDfsdT << " 1/K" << nl
            << "    CpApp min/max         = "
            << minCpApp << " " << maxCpApp << " J/(kg K)" << nl
            << "    liquid/mushy/solid    = "
            << nLiquid << " " << nMushy << " " << nSolid << nl
            << "    max mass-closure err  = "
            << maxMassClosureError << nl
            << "    max |delta fs|        = "
            << maxDeltaFs
            << endl;
    }

    return maxDeltaFs;
}


void Foam::solvers::continuousCastingMacrosegregation::updateBuoyancy()
{
    const dimensionedVector gDim("g", g_.dimensions(), g_.value());

    const dimensionedScalar betaTDim
    (
        "betaT",
        dimless/dimTemperature,
        betaT_
    );

    const dimensionedScalar betaCDim
    (
        "betaC",
        dimless,
        betaC_
    );

    const dimensionedScalar TRefDim
    (
        "TRef",
        dimTemperature,
        TRef_
    );

    const dimensionedScalar C0Dim
    (
        "Carbon0",
        dimless,
        Carbon0_
    );

    // Dynamic Boussinesq acceleration after absorbing the hydrostatic
    // reference contribution into pressure:
    //
    // a_b,T = -g betaT (T - T0)
    // a_b,C = -g betaC (Cl - C0)
    //
    // No phase-fraction weighting is used here. This follows the benchmark
    // momentum model; BKC drag will separately suppress motion in mush/solid.

    buoyancyThermal_ =
        -gDim*betaTDim*(T_ - TRefDim);

    buoyancySolutal_ =
        -gDim*betaCDim*(CarbonL_ - C0Dim);

    buoyancyTotal_ =
        buoyancyThermal_ + buoyancySolutal_;

    buoyancyThermal_.correctBoundaryConditions();
    buoyancySolutal_.correctBoundaryConditions();
    buoyancyTotal_.correctBoundaryConditions();

    scalar minDeltaT = GREAT;
    scalar maxDeltaT = -GREAT;
    scalar minDeltaC = GREAT;
    scalar maxDeltaC = -GREAT;

    scalar maxAThermal = 0.0;
    scalar maxASolutal = 0.0;
    scalar maxATotal = 0.0;

    scalar integralFThermal = 0.0;
    scalar integralFSolutal = 0.0;

    vector netFThermal = vector::zero;
    vector netFSolutal = vector::zero;
    vector netFTotal = vector::zero;

    const scalarField& V = mesh_.V();

    forAll(T_, celli)
    {
        const scalar dT = T_[celli] - TRef_;
        const scalar dC = CarbonL_[celli] - Carbon0_;

        minDeltaT = min(minDeltaT, dT);
        maxDeltaT = max(maxDeltaT, dT);
        minDeltaC = min(minDeltaC, dC);
        maxDeltaC = max(maxDeltaC, dC);

        const vector aT = buoyancyThermal_[celli];
        const vector aC = buoyancySolutal_[celli];
        const vector aTot = buoyancyTotal_[celli];

        maxAThermal = max(maxAThermal, mag(aT));
        maxASolutal = max(maxASolutal, mag(aC));
        maxATotal = max(maxATotal, mag(aTot));

        integralFThermal += rho_*mag(aT)*V[celli];
        integralFSolutal += rho_*mag(aC)*V[celli];

        netFThermal += rho_*aT*V[celli];
        netFSolutal += rho_*aC*V[celli];
        netFTotal += rho_*aTot*V[celli];
    }

    minDeltaT = returnReduce(minDeltaT, minOp<scalar>());
    maxDeltaT = returnReduce(maxDeltaT, maxOp<scalar>());
    minDeltaC = returnReduce(minDeltaC, minOp<scalar>());
    maxDeltaC = returnReduce(maxDeltaC, maxOp<scalar>());

    maxAThermal = returnReduce(maxAThermal, maxOp<scalar>());
    maxASolutal = returnReduce(maxASolutal, maxOp<scalar>());
    maxATotal = returnReduce(maxATotal, maxOp<scalar>());

    integralFThermal =
        returnReduce(integralFThermal, sumOp<scalar>());

    integralFSolutal =
        returnReduce(integralFSolutal, sumOp<scalar>());

    netFThermal =
        returnReduce(netFThermal, sumOp<vector>());

    netFSolutal =
        returnReduce(netFSolutal, sumOp<vector>());

    netFTotal =
        returnReduce(netFTotal, sumOp<vector>());

    Info<< "Thermo-solutal buoyancy audit" << nl
        << "    deltaT min/max          = "
        << minDeltaT << " " << maxDeltaT << " K" << nl
        << "    deltaCarbonL min/max    = "
        << minDeltaC << " " << maxDeltaC << nl
        << "    max(|a_b,T|)            = "
        << maxAThermal << " m/s2" << nl
        << "    max(|a_b,C|)            = "
        << maxASolutal << " m/s2" << nl
        << "    max(|a_b,total|)        = "
        << maxATotal << " m/s2" << nl
        << "    integral(|F_b,T| dV)   = "
        << integralFThermal << " N" << nl
        << "    integral(|F_b,C| dV)   = "
        << integralFSolutal << " N" << nl
        << "    net thermal force       = "
        << netFThermal << " N" << nl
        << "    net solutal force       = "
        << netFSolutal << " N" << nl
        << "    net total force         = "
        << netFTotal << " N"
        << endl;
}


void Foam::solvers::continuousCastingMacrosegregation::updateBKCDrag
(
    const bool report
)
{
    // Standard BKC, paper Eqs. (5) and (7):
    //
    //   K^-1 = K0^-1 fsB^2/(1-fsB)^3
    //   K0   = lambda2^2/180
    //   fsB  = min(fs, 0.99)
    //
    // us = 0 for the AFRODITE benchmark. After division of momentum by
    // constant rho, the Darcy acceleration is
    //
    //   a_D = -(mu_l/rho) K^-1 U = -(nu_l/K) U.
    //
    // The standard-BKC model omits the relative-advection momentum term Sr.

    const scalar K0 = sqr(lambda2_)/180.0;
    const scalar nuLiquid = muLiquid_/rho_;

    scalar minInvK = GREAT;
    scalar maxInvK = -GREAT;
    scalar minDragCoeff = GREAT;
    scalar maxDragCoeff = -GREAT;
    scalar maxDragAcceleration = 0.0;
    scalar maxUMushy = 0.0;
    scalar maxUSolidLike = 0.0;
    scalar integralDragForce = 0.0;
    vector netDragForce = vector::zero;

    const scalarField& V = mesh_.V();

    forAll(fs_, celli)
    {
        const scalar fsCell =
            min(max(fs_[celli], scalar(0)), scalar(1));

        const scalar fsB =
            min(fsCell, scalar(0.99));

        scalar invK = 0.0;

        if (fsB > SMALL)
        {
            invK =
                sqr(fsB)
               /(K0*pow3(1.0 - fsB));
        }

        const scalar dragCoeff =
            nuLiquid*invK;

        const vector dragAcceleration =
            -dragCoeff*U_[celli];

        bkcInvK_[celli] = invK;
        bkcDragCoeff_[celli] = dragCoeff;
        bkcDragAcceleration_[celli] = dragAcceleration;

        minInvK = min(minInvK, invK);
        maxInvK = max(maxInvK, invK);
        minDragCoeff = min(minDragCoeff, dragCoeff);
        maxDragCoeff = max(maxDragCoeff, dragCoeff);

        maxDragAcceleration =
            max(maxDragAcceleration, mag(dragAcceleration));

        if (fsCell > SMALL && fsCell < 1.0 - SMALL)
        {
            maxUMushy =
                max(maxUMushy, mag(U_[celli]));
        }

        if (fsCell >= 0.99)
        {
            maxUSolidLike =
                max(maxUSolidLike, mag(U_[celli]));
        }

        const vector dragForce =
            rho_*dragAcceleration;

        integralDragForce +=
            mag(dragForce)*V[celli];

        netDragForce +=
            dragForce*V[celli];
    }

    bkcInvK_.correctBoundaryConditions();
    bkcDragCoeff_.correctBoundaryConditions();
    bkcDragAcceleration_.correctBoundaryConditions();

    if (!report)
    {
        return;
    }

    minInvK = returnReduce(minInvK, minOp<scalar>());
    maxInvK = returnReduce(maxInvK, maxOp<scalar>());
    minDragCoeff = returnReduce(minDragCoeff, minOp<scalar>());
    maxDragCoeff = returnReduce(maxDragCoeff, maxOp<scalar>());
    maxDragAcceleration =
        returnReduce(maxDragAcceleration, maxOp<scalar>());
    maxUMushy = returnReduce(maxUMushy, maxOp<scalar>());
    maxUSolidLike = returnReduce(maxUSolidLike, maxOp<scalar>());
    integralDragForce =
        returnReduce(integralDragForce, sumOp<scalar>());
    netDragForce =
        returnReduce(netDragForce, sumOp<vector>());

    Info<< "Standard-BKC permeability audit" << nl
        << "    K0                     = "
        << K0 << " m2" << nl
        << "    fsB cap                = 0.99" << nl
        << "    invK min/max           = "
        << minInvK << " " << maxInvK << " 1/m2" << nl
        << "    dragCoeff min/max      = "
        << minDragCoeff << " " << maxDragCoeff << " 1/s" << nl
        << "    max(|a_D|)             = "
        << maxDragAcceleration << " m/s2" << nl
        << "    max(|U|) mushy         = "
        << maxUMushy << " m/s" << nl
        << "    max(|U|) fs>=0.99      = "
        << maxUSolidLike << " m/s" << nl
        << "    integral(|F_D| dV)     = "
        << integralDragForce << " N" << nl
        << "    net drag force         = "
        << netDragForce << " N"
        << endl;
}


void Foam::solvers::continuousCastingMacrosegregation::updateSpeciesDiagnostics
(
    const bool report
)
{
    scalar minD = GREAT;
    scalar maxD = -GREAT;
    scalar minQ = GREAT;
    scalar maxQ = -GREAT;
    scalar minMS = GREAT;
    scalar maxMS = -GREAT;

    scalar minC = GREAT;
    scalar maxC = -GREAT;
    scalar minCl = GREAT;
    scalar maxCl = -GREAT;
    scalar minCs = GREAT;
    scalar maxCs = -GREAT;

    scalar domainVolume = 0.0;
    scalar soluteInventory = 0.0;

    const scalarField& V = mesh_.V();

    forAll(Carbon_, celli)
    {
        const scalar fsCell =
            min(max(fs_[celli], scalar(0)), scalar(1));

        const scalar flCell =
            1.0 - fsCell;

        const scalar DCell =
            fsCell*DS_ + flCell*DL_;

        const scalar denominator =
            1.0 + fsCell*(kp_ - 1.0);

        if (denominator <= SMALL)
        {
            FatalErrorInFunction
                << "Invalid Lever-rule denominator while assembling species "
                << "transport." << nl
                << "    cell = " << celli << nl
                << "    fs   = " << fsCell << nl
                << "    kp   = " << kp_
                << abort(FatalError);
        }

        const scalar qCell =
            1.0/denominator;

        const scalar MSCell =
            100.0*(Carbon_[celli] - Carbon0_)/Carbon0_;

        speciesDiffusivity_[celli] = DCell;
        liquidAdvectionFactor_[celli] = qCell;
        macrosegregation_[celli] = MSCell;

        minD = min(minD, DCell);
        maxD = max(maxD, DCell);

        minQ = min(minQ, qCell);
        maxQ = max(maxQ, qCell);

        minMS = min(minMS, MSCell);
        maxMS = max(maxMS, MSCell);

        minC = min(minC, Carbon_[celli]);
        maxC = max(maxC, Carbon_[celli]);

        minCl = min(minCl, CarbonL_[celli]);
        maxCl = max(maxCl, CarbonL_[celli]);

        minCs = min(minCs, CarbonS_[celli]);
        maxCs = max(maxCs, CarbonS_[celli]);

        domainVolume += V[celli];
        soluteInventory += Carbon_[celli]*V[celli];
    }

    speciesDiffusivity_.correctBoundaryConditions();
    liquidAdvectionFactor_.correctBoundaryConditions();
    macrosegregation_.correctBoundaryConditions();

    if (!report)
    {
        return;
    }

    minD = returnReduce(minD, minOp<scalar>());
    maxD = returnReduce(maxD, maxOp<scalar>());

    minQ = returnReduce(minQ, minOp<scalar>());
    maxQ = returnReduce(maxQ, maxOp<scalar>());

    minMS = returnReduce(minMS, minOp<scalar>());
    maxMS = returnReduce(maxMS, maxOp<scalar>());

    minC = returnReduce(minC, minOp<scalar>());
    maxC = returnReduce(maxC, maxOp<scalar>());

    minCl = returnReduce(minCl, minOp<scalar>());
    maxCl = returnReduce(maxCl, maxOp<scalar>());

    minCs = returnReduce(minCs, minOp<scalar>());
    maxCs = returnReduce(maxCs, maxOp<scalar>());

    domainVolume =
        returnReduce(domainVolume, sumOp<scalar>());

    soluteInventory =
        returnReduce(soluteInventory, sumOp<scalar>());

    const scalar meanCarbon =
        soluteInventory/domainVolume;

    const scalar referenceInventory =
        Carbon0_*domainVolume;

    const scalar relativeInventoryError =
        (soluteInventory - referenceInventory)
       /max(mag(referenceInventory), SMALL);

    Info<< "Mixture-species transport audit" << nl
        << "    Carbon min/max          = "
        << minC << " " << maxC << nl
        << "    CarbonL min/max         = "
        << minCl << " " << maxCl << nl
        << "    CarbonS min/max         = "
        << minCs << " " << maxCs << nl
        << "    Dmix min/max            = "
        << minD << " " << maxD << " m2/s" << nl
        << "    CarbonL/Carbon min/max  = "
        << minQ << " " << maxQ << nl
        << "    MS min/max              = "
        << minMS << " " << maxMS << " %" << nl
        << "    volumeAverage(Carbon)   = "
        << meanCarbon << nl
        << "    integral(Carbon dV)     = "
        << soluteInventory << " m3" << nl
        << "    relative inventory err  = "
        << relativeInventoryError
        << endl;
}


Foam::scalar
Foam::solvers::continuousCastingMacrosegregation::solveSpeciesTransport()
{
    // Bennon-Incropera mixture species equation used by the paper, Eq. (19):
    //
    // d(C)/dt + div(U C)
    //   = div(D grad(C))
    //   + div(fl Dl grad(Cl-C))
    //   + div(fs Ds grad(Cs-C))
    //   - div((U-us)(Cl-C)).
    //
    // For the AFRODITE standard-BKC model us = 0.
    //
    // IMPORTANT v9 CHANGE
    // -------------------
    // v8 analytically combined the two advective terms into div(U Cl) and
    // discretized that collapsed flux implicitly using Cl = q*C.
    //
    // Although the continuous equations are algebraically equivalent, the
    // finite-volume operators need not be numerically identical when using
    // bounded linear-upwind interpolation. v9 therefore discretizes Eq. (19)
    // term-by-term, following the published equation directly:
    //
    //   + div(U C)              : implicit in Carbon
    //   - div(U (Cl-C))         : explicit RHS correction
    //
    // The mixture diffusion term is implicit. The liquid/solid
    // phase-difference diffusion terms remain explicit. For the published
    // AFRODITE properties DS = 0.

    updateSpeciesDiagnostics(false);

    const dimensionedScalar DLDim
    (
        "DL",
        dimArea/dimTime,
        DL_
    );

    const dimensionedScalar DSDim
    (
        "DS",
        dimArea/dimTime,
        DS_
    );

    const volScalarField liquidFraction
    (
        IOobject
        (
            "liquidFractionSpeciesTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        scalar(1) - fs_
    );

    // Explicit relative-composition field appearing in the final advective
    // term of Eq. (19). For standard BKC us=0.
    const volScalarField relativeLiquidComposition
    (
        IOobject
        (
            "relativeLiquidCompositionTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        CarbonL_ - Carbon_
    );

    const scalarField CarbonBefore(Carbon_.primitiveField());

    fvScalarMatrix CarbonEqn
    (
        fvm::ddt(Carbon_)
      + fvm::div(phi_, Carbon_, "div(phi,Carbon)")
      - fvm::laplacian(speciesDiffusivity_, Carbon_)
     ==
        fvc::laplacian
        (
            liquidFraction*DLDim,
            CarbonL_ - Carbon_
        )
      + fvc::laplacian
        (
            fs_*DSDim,
            CarbonS_ - Carbon_
        )
      - fvc::div
        (
            phi_,
            relativeLiquidComposition,
            "div(phi,Carbon)"
        )
    );

    CarbonEqn.solve();

    scalar maxDeltaCarbon = 0.0;

    forAll(Carbon_, celli)
    {
        maxDeltaCarbon =
            max
            (
                maxDeltaCarbon,
                mag(Carbon_[celli] - CarbonBefore[celli])
            );
    }

    maxDeltaCarbon =
        returnReduce(maxDeltaCarbon, maxOp<scalar>());

    return maxDeltaCarbon;
}


void Foam::solvers::continuousCastingMacrosegregation::updateEnergyDiagnostics() const
{
    // Eq. (13), standard BKC, us = 0:
    //
    //   d(rho Cp T)/dt + div(rho U Cp T)
    //     = div(k grad T) + d(rho L fs)/dt.
    //
    // In the closed AFRODITE cavity the normal advective flux is zero.
    // Rearranging the latent term gives the globally conserved quantity
    //
    //   E* = integral rho (Cp*T - L*fs) dV,
    //
    // for which
    //
    //   dE*/dt = integral_boundary k grad(T).n dA.
    //
    // The sign convention below is therefore positive when the boundary
    // conduction term adds energy to the domain and negative when heat leaves.

    const scalar deltaTValue = runTime.deltaTValue();

    if (deltaTValue <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive time step in energy audit: "
            << deltaTValue
            << abort(FatalError);
    }

    scalar sensibleEnergy = 0.0;
    scalar latentReferenceEnergy = 0.0;
    scalar effectiveEnergy = 0.0;

    scalar oldSensibleEnergy = 0.0;
    scalar oldLatentReferenceEnergy = 0.0;
    scalar oldEffectiveEnergy = 0.0;

    const scalarField& V = mesh_.V();

    const volScalarField& CpOld = CpMix_.oldTime();
    const volScalarField& TOld = T_.oldTime();
    const volScalarField& fsOld = fs_.oldTime();

    forAll(T_, celli)
    {
        const scalar sensibleCell =
            rho_*CpMix_[celli]*T_[celli]*V[celli];

        const scalar latentCell =
            rho_*latentHeat_*fs_[celli]*V[celli];

        const scalar oldSensibleCell =
            rho_*CpOld[celli]*TOld[celli]*V[celli];

        const scalar oldLatentCell =
            rho_*latentHeat_*fsOld[celli]*V[celli];

        sensibleEnergy += sensibleCell;
        latentReferenceEnergy += latentCell;
        effectiveEnergy += sensibleCell - latentCell;

        oldSensibleEnergy += oldSensibleCell;
        oldLatentReferenceEnergy += oldLatentCell;
        oldEffectiveEnergy += oldSensibleCell - oldLatentCell;
    }

    sensibleEnergy =
        returnReduce(sensibleEnergy, sumOp<scalar>());

    latentReferenceEnergy =
        returnReduce(latentReferenceEnergy, sumOp<scalar>());

    effectiveEnergy =
        returnReduce(effectiveEnergy, sumOp<scalar>());

    oldSensibleEnergy =
        returnReduce(oldSensibleEnergy, sumOp<scalar>());

    oldLatentReferenceEnergy =
        returnReduce(oldLatentReferenceEnergy, sumOp<scalar>());

    oldEffectiveEnergy =
        returnReduce(oldEffectiveEnergy, sumOp<scalar>());

    scalar conductiveBoundaryPower = 0.0;

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        // Processor/cyclic-type coupled patches are internal to the global
        // domain and must not contribute to the physical boundary balance.
        if (patch.coupled())
        {
            continue;
        }

        const fvPatchScalarField& Tp =
            T_.boundaryField()[patchi];

        const fvPatchScalarField& kp =
            kEff_.boundaryField()[patchi];

        const scalarField& magSf = patch.magSf();

        const scalarField snGradT(Tp.snGrad());

        forAll(snGradT, facei)
        {
            conductiveBoundaryPower +=
                kp[facei]*snGradT[facei]*magSf[facei];
        }
    }

    conductiveBoundaryPower =
        returnReduce(conductiveBoundaryPower, sumOp<scalar>());

    const scalar effectiveEnergyRate =
        (effectiveEnergy - oldEffectiveEnergy)/deltaTValue;

    const scalar balanceResidual =
        effectiveEnergyRate - conductiveBoundaryPower;

    const scalar balanceScale =
        max
        (
            max
            (
                mag(effectiveEnergyRate),
                mag(conductiveBoundaryPower)
            ),
            SMALL
        );

    const scalar relativeBalanceResidual =
        balanceResidual/balanceScale;

    Info<< "Global energy-conservation audit" << nl
        << "    sensible energy             = "
        << sensibleEnergy << " J" << nl
        << "    latent reference rho*L*fs  = "
        << latentReferenceEnergy << " J" << nl
        << "    effective E*=Esens-Elat    = "
        << effectiveEnergy << " J" << nl
        << "    old effective E*           = "
        << oldEffectiveEnergy << " J" << nl
        << "    dE*/dt                     = "
        << effectiveEnergyRate << " W" << nl
        << "    boundary integral k*dT/dn  = "
        << conductiveBoundaryPower << " W" << nl
        << "    energy balance residual    = "
        << balanceResidual << " W" << nl
        << "    relative balance residual  = "
        << relativeBalanceResidual
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * //

Foam::solvers::continuousCastingMacrosegregation::continuousCastingMacrosegregation
(
    fvMesh& mesh
)
:
    basicFluidSolver(mesh),

    p_
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensions::kinematicPressure
    ),

    pressureReference(p_, pimple.dict()),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensions::velocity
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & mesh.Sf()
    ),

    alloyProperties_
    (
        IOobject
        (
            "alloyProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    rho_(alloyProperties_.lookup<scalar>("rho")),
    CpLiquid_(alloyProperties_.lookup<scalar>("CpLiquid")),
    CpSolid_(alloyProperties_.lookup<scalar>("CpSolid")),
    kLiquid_(alloyProperties_.lookup<scalar>("kLiquid")),
    kSolid_(alloyProperties_.lookup<scalar>("kSolid")),
    latentHeat_(alloyProperties_.lookup<scalar>("latentHeat")),
    Carbon0_(alloyProperties_.lookup<scalar>("Carbon0")),
    kp_(alloyProperties_.lookup<scalar>("kp")),
    liquidusSlope_(alloyProperties_.lookup<scalar>("liquidusSlope")),
    Tmelt_(alloyProperties_.lookup<scalar>("Tmelt")),
    muLiquid_(alloyProperties_.lookup<scalar>("muLiquid")),
    DL_(alloyProperties_.lookup<scalar>("DL")),
    DS_(alloyProperties_.lookup<scalar>("DS")),
    betaT_(alloyProperties_.lookup<scalar>("betaT")),
    betaC_(alloyProperties_.lookup<scalar>("betaC")),
    TRef_(alloyProperties_.lookup<scalar>("TRef")),
    lambda2_(alloyProperties_.lookup<scalar>("lambda2")),
    solidVelocity_
    (
        alloyProperties_.lookupOrDefault<vector>
        (
            "solidVelocity",
            vector::zero
        )
    ),
    nSolidificationLoops_
    (
        alloyProperties_.lookup<label>("nSolidificationLoops")
    ),

    g_
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    T_
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Carbon_
    (
        IOobject
        (
            "Carbon",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    fs_
    (
        IOobject
        (
            "fs",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    CarbonL_
    (
        IOobject
        (
            "CarbonL",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    CarbonS_
    (
        IOobject
        (
            "CarbonS",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    speciesDiffusivity_
    (
        IOobject
        (
            "speciesDiffusivity",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DL", dimArea/dimTime, DL_),
        zeroGradientFvPatchScalarField::typeName
    ),

    liquidAdvectionFactor_
    (
        IOobject
        (
            "liquidAdvectionFactor",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    macrosegregation_
    (
        IOobject
        (
            "macrosegregation",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    CpMix_
    (
        IOobject
        (
            "CpMix",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "CpLiquid",
            dimensionSet(0, 2, -2, -1, 0, 0, 0),
            CpLiquid_
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    kEff_
    (
        IOobject
        (
            "kEff",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "kLiquid",
            dimensionSet(1, 1, -3, -1, 0, 0, 0),
            kLiquid_
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    dfsdT_
    (
        IOobject
        (
            "dfsdT",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimensionSet(0, 0, 0, -1, 0, 0, 0),
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    CpApp_
    (
        IOobject
        (
            "CpApp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Cp",
            dimensionSet(0, 2, -2, -1, 0, 0, 0),
            CpLiquid_
        ),
        zeroGradientFvPatchScalarField::typeName
    ),

    buoyancyThermal_
    (
        IOobject
        (
            "buoyancyThermal",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength/sqr(dimTime), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    buoyancySolutal_
    (
        IOobject
        (
            "buoyancySolutal",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength/sqr(dimTime), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    buoyancyTotal_
    (
        IOobject
        (
            "buoyancyTotal",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength/sqr(dimTime), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    bkcInvK_
    (
        IOobject
        (
            "bkcInvK",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(sqr(dimLength)), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    bkcDragCoeff_
    (
        IOobject
        (
            "bkcDragCoeff",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    bkcDragAcceleration_
    (
        IOobject
        (
            "bkcDragAcceleration",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength/sqr(dimTime), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    ),

    viscosity(viscosityModel::New(mesh)),

    momentumTransport
    (
        incompressible::momentumTransportModel::New
        (
            U_,
            phi_,
            viscosity
        )
    ),

    MRF(MRFZones::New(mesh)),

    p(p_),
    U(U_),
    phi(phi_),
    T(T_),
    Carbon(Carbon_),
    fs(fs_),
    CarbonL(CarbonL_),
    CarbonS(CarbonS_),
    speciesDiffusivity(speciesDiffusivity_),
    liquidAdvectionFactor(liquidAdvectionFactor_),
    macrosegregation(macrosegregation_),
    CpMix(CpMix_),
    kEff(kEff_),
    dfsdT(dfsdT_),
    CpApp(CpApp_),
    bkcInvK(bkcInvK_),
    bkcDragCoeff(bkcDragCoeff_),
    bkcDragAcceleration(bkcDragAcceleration_)
{
    validateAlloyProperties();
    updatePhaseState(true);
    updateSpeciesDiagnostics(true);

    const scalar alphaTLiquid =
        kLiquid_/(rho_*CpLiquid_);

    const scalar alphaTSolid =
        kSolid_/(rho_*CpSolid_);

    Info<< "continuousCastingMacrosegregation: species discretization" << nl
        << "    Eq. (19) advection = term-by-term" << nl
        << "    div(U*C)           = implicit" << nl
        << "    div(U*(Cl-C))      = explicit RHS correction"
        << endl;

    Info<< "continuousCastingMacrosegregation: energy discretization" << nl
        << "    Eq. (13) advection = direct CpMix*T product" << nl
        << "    implementation     = implicit split + deferred correction" << nl
        << "    v11c change        = diagnostics only; solved equations unchanged"
        << endl;

    Info<< "continuousCastingMacrosegregation: solid-velocity infrastructure" << nl
        << "    solidVelocity     = " << solidVelocity_ << " m/s" << nl
        << "    |solidVelocity|   = " << mag(solidVelocity_) << " m/s" << nl
        << "    CC-1 coupling     = diagnostics only; equations unchanged"
        << endl;

    Info<< "continuousCastingMacrosegregation: alloy properties" << nl
        << "    rho             = " << rho_ << " kg/m3" << nl
        << "    CpLiquid        = " << CpLiquid_ << " J/(kg K)" << nl
        << "    CpSolid         = " << CpSolid_ << " J/(kg K)" << nl
        << "    kLiquid         = " << kLiquid_ << " W/(m K)" << nl
        << "    kSolid          = " << kSolid_ << " W/(m K)" << nl
        << "    alphaTLiquid    = " << alphaTLiquid << " m2/s" << nl
        << "    alphaTSolid     = " << alphaTSolid << " m2/s" << nl
        << "    latentHeat      = " << latentHeat_ << " J/kg" << nl
        << "    Carbon0         = " << Carbon0_ << nl
        << "    kp              = " << kp_ << nl
        << "    liquidusSlope   = " << liquidusSlope_
        << " K/(mass fraction)" << nl
        << "    Tmelt           = " << Tmelt_ << " K" << nl
        << "    muLiquid        = " << muLiquid_ << " Pa s" << nl
        << "    DL              = " << DL_ << " m2/s" << nl
        << "    DS              = " << DS_ << " m2/s" << nl
        << "    betaT           = " << betaT_ << " 1/K" << nl
        << "    betaC           = " << betaC_
        << " 1/(solute mass fraction)" << nl
        << "    TRef            = " << TRef_ << " K" << nl
        << "    g               = " << g_.value() << " m/s2" << nl
        << "    lambda2         = " << lambda2_ << " m" << nl
        << "    solidification loops = "
        << nSolidificationLoops_
        << endl;

    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * //

Foam::solvers::continuousCastingMacrosegregation::~continuousCastingMacrosegregation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::continuousCastingMacrosegregation::preSolve()
{
    if ((mesh.dynamic() || MRF.size()) && !Uf.valid())
    {
        Info<< "Constructing face momentum Uf" << endl;

        // Ensure the U BCs are up-to-date before constructing Uf
        U_.correctBoundaryConditions();

        Uf = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U)
        );
    }

    fvModels().preUpdateMesh();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::continuousCastingMacrosegregation::prePredictor()
{}


void Foam::solvers::continuousCastingMacrosegregation::momentumTransportPredictor()
{
    momentumTransport->predict();
}


void Foam::solvers::continuousCastingMacrosegregation::
thermophysicalTransportPredictor()
{}


void Foam::solvers::continuousCastingMacrosegregation::thermophysicalPredictor()
{
    // Standard-BKC energy equation following Moeinirad & Amani Eq. (13).
    //
    // For standard BKC, us = 0, so the two convective latent-enthalpy terms
    // cancel exactly:
    //
    //   +div(rho L fs U) - div(rho L fs (U-us)) = 0.
    //
    // With constant rho, but phase-dependent mixture Cp and k:
    //
    //   d(CpMix*T)/dt + div(U CpMix T)
    //       = div((kEff/rho) grad(T)) + L d(fs)/dt.
    //
    // Eq. (16) recasts the local latent term as
    //
    //   d(fs)/dt = (dfs/dT) dT/dt,
    //
    // with dfs/dT evaluated analytically from the Lever rule. Therefore:
    //
    //   d(CpMix*T)/dt
    // + [-L dfs/dT] dT/dt
    // + div(U CpMix T)
    // - div((kEff/rho) grad(T)) = 0.
    //
    // CpMix and kEff are frozen within each fixed-point solidification
    // correction and are updated from the new phase state after each solve.

    const dimensionedScalar rhoDim
    (
        "rho",
        dimDensity,
        rho_
    );

    const dimensionedScalar latentHeatDim
    (
        "latentHeat",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        latentHeat_
    );

    const scalar deltaTValue = runTime.deltaTValue();

    if (deltaTValue <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive time step: " << deltaTValue
            << abort(FatalError);
    }

    const dimensionedScalar rDeltaT
    (
        "rDeltaT",
        dimless/dimTime,
        1.0/deltaTValue
    );

    // Preserve previous physical-time algebraic fields while their current
    // values are changed during the inner solidification iterations.
    // These old-time fields are also required by the global Eq. (13)
    // energy-conservation audit.
    CpMix_.oldTime();
    fs_.oldTime();

    Info<< "Solidification/variable-property energy coupling" << nl
        << "    Eq. (13) sensible advection = direct CpMix*T product" << nl
        << "    loops = " << nSolidificationLoops_
        << endl;

    for (label corr = 0; corr < nSolidificationLoops_; ++corr)
    {
        // Positive in the mush because dfs/dT < 0.
        const tmp<volScalarField> tLatentCp =
            -latentHeatDim*dfsdT_;

        // Effective conductivity term after dividing Eq. (13) by rho.
        const tmp<volScalarField> tKbyRho =
            kEff_/rhoDim;

        // Eq. (13) contains div(U * CpMix * T).
        //
        // v9 assembled this as interpolate(CpMix)*phi acting implicitly on T.
        // That is algebraically equivalent in the continuous equation, but
        // not discretely identical to applying bounded linear-upwind to the
        // transported product CpMix*T.
        //
        // v10 keeps the stable implicit split operator as a matrix base and
        // applies a deferred explicit correction:
        //
        //   directAdv = div(phi, CpMix*T)
        //   splitAdv  = div(interpolate(CpMix)*phi, T)
        //   correction = directAdv - splitAdv
        //
        // The correction is placed on the RHS with a minus sign, so at fixed-
        // point convergence the effective discrete advection operator is the
        // direct product operator appearing in Eq. (13).

        const volScalarField CpT
        (
            IOobject
            (
                "CpT",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            CpMix_*T_
        );

        const surfaceScalarField phiCp
        (
            IOobject
            (
                "phiCp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvc::interpolate(CpMix_)*phi_
        );

        const tmp<volScalarField> tDirectEnergyAdvection =
            fvc::div(phi_, CpT, "div(phi,T)");

        const tmp<volScalarField> tSplitEnergyAdvection =
            fvc::div(phiCp, T_, "div(phi,T)");

        const tmp<volScalarField> tEnergyAdvectionCorrection =
            tDirectEnergyAdvection - tSplitEnergyAdvection;

        scalar maxEnergyAdvectionCorrection = 0.0;

        const volScalarField& energyAdvectionCorrection =
            tEnergyAdvectionCorrection();

        forAll(energyAdvectionCorrection, celli)
        {
            maxEnergyAdvectionCorrection =
                max
                (
                    maxEnergyAdvectionCorrection,
                    mag(energyAdvectionCorrection[celli])
                );
        }

        maxEnergyAdvectionCorrection =
            returnReduce
            (
                maxEnergyAdvectionCorrection,
                maxOp<scalar>()
            );

        const scalarField TBefore(T_.primitiveField());

        fvScalarMatrix TEqn
        (
            fvm::ddt(CpMix_, T_)
          + fvm::Sp(rDeltaT*tLatentCp(), T_)
          + fvm::div(phiCp, T_, "div(phi,T)")
          - fvm::laplacian(tKbyRho(), T_)
         ==
            rDeltaT*tLatentCp()*T_.oldTime()
          - tEnergyAdvectionCorrection()
        );

        TEqn.solve();

        // -----------------------------------------------------------------
        // v11 discrete Eq. (13) audit -- DIAGNOSTIC ONLY
        //
        // This audit is evaluated immediately after the T solve and before
        // Carbon/phase-equilibrium updates change CpMix, kEff or dfsdT.
        // Therefore it uses the same frozen coefficients and the same
        // deferred-correction term that actually entered TEqn.
        //
        // The solved equation, after division by rho, is:
        //
        //   d(Cp*T)/dt
        // + (-L*dfs/dT)*(T-Told)/dt
        // + div(interpolate(Cp)*phi,T)
        // + [directAdv(old iterate) - splitAdv(old iterate)]
        // - laplacian(k/rho,T)
        // = 0.
        //
        // Multiplication by rho*V and global summation gives watts.
        //
        // No field used by the solver is modified below.

        const bool finalThermalCorrection =
            corr == nSolidificationLoops_ - 1;

        if (finalThermalCorrection)
        {
            // Foundation v14 has no fvc::ddt(CpMix_, T_) overload.
            // For this fixed mesh with Euler time integration, the explicit
            // diagnostic equivalent of fvm::ddt(CpMix_, T_) is
            //
            //   rDeltaT*(Cp*T - Cp.oldTime()*T.oldTime()).
            //
            // This is diagnostic-only and does not alter TEqn.
            const tmp<volScalarField> tDiscreteSensibleStorage =
                rDeltaT
               *(
                    CpMix_*T_
                  - CpMix_.oldTime()*T_.oldTime()
                );

            const tmp<volScalarField> tDiscreteLatentStorage =
                rDeltaT*tLatentCp()*(T_ - T_.oldTime());

            const tmp<volScalarField> tDiscreteSplitAdvection =
                fvc::div(phiCp, T_, "div(phi,T)");

            // This is the correction actually used on the RHS of TEqn.
            const volScalarField& correctionUsed =
                tEnergyAdvectionCorrection();

            const tmp<volScalarField> tDiscreteEffectiveAdvection =
                tDiscreteSplitAdvection() + correctionUsed;

            const tmp<volScalarField> tDiscreteConduction =
                fvc::laplacian(tKbyRho(), T_);

            const tmp<volScalarField> tDiscreteResidual =
                tDiscreteSensibleStorage()
              + tDiscreteLatentStorage()
              + tDiscreteEffectiveAdvection()
              - tDiscreteConduction();

            // Also evaluate the direct Cp*T advection with the converged T.
            // Its difference from the operator actually used measures the
            // remaining lag introduced by the deferred correction.
            const volScalarField CpTFinal
            (
                IOobject
                (
                    "CpTFinalAudit",
                    runTime.name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                CpMix_*T_
            );

            const tmp<volScalarField> tDirectAdvectionFinal =
                fvc::div(phi_, CpTFinal, "div(phi,T)");

            const tmp<volScalarField> tDeferredLag =
                tDirectAdvectionFinal()
              - tDiscreteEffectiveAdvection();

            scalar sensibleStoragePower = 0.0;
            scalar latentStoragePower = 0.0;
            scalar effectiveAdvectionPower = 0.0;
            scalar conductionPower = 0.0;
            scalar discreteResidualPower = 0.0;
            scalar directAdvectionFinalPower = 0.0;
            scalar deferredLagPower = 0.0;

            scalar maxDiscreteResidualPerMass = 0.0;
            scalar maxDeferredLagPerMass = 0.0;

            const scalarField& V = mesh_.V();

            const volScalarField& sensibleStorage =
                tDiscreteSensibleStorage();

            const volScalarField& latentStorage =
                tDiscreteLatentStorage();

            const volScalarField& effectiveAdvection =
                tDiscreteEffectiveAdvection();

            const volScalarField& conduction =
                tDiscreteConduction();

            const volScalarField& discreteResidual =
                tDiscreteResidual();

            const volScalarField& directAdvectionFinal =
                tDirectAdvectionFinal();

            const volScalarField& deferredLag =
                tDeferredLag();

            forAll(T_, celli)
            {
                const scalar rhoV = rho_*V[celli];

                sensibleStoragePower +=
                    rhoV*sensibleStorage[celli];

                latentStoragePower +=
                    rhoV*latentStorage[celli];

                effectiveAdvectionPower +=
                    rhoV*effectiveAdvection[celli];

                conductionPower +=
                    rhoV*conduction[celli];

                discreteResidualPower +=
                    rhoV*discreteResidual[celli];

                directAdvectionFinalPower +=
                    rhoV*directAdvectionFinal[celli];

                deferredLagPower +=
                    rhoV*deferredLag[celli];

                maxDiscreteResidualPerMass =
                    max
                    (
                        maxDiscreteResidualPerMass,
                        mag(discreteResidual[celli])
                    );

                maxDeferredLagPerMass =
                    max
                    (
                        maxDeferredLagPerMass,
                        mag(deferredLag[celli])
                    );
            }

            sensibleStoragePower =
                returnReduce(sensibleStoragePower, sumOp<scalar>());

            latentStoragePower =
                returnReduce(latentStoragePower, sumOp<scalar>());

            effectiveAdvectionPower =
                returnReduce(effectiveAdvectionPower, sumOp<scalar>());

            conductionPower =
                returnReduce(conductionPower, sumOp<scalar>());

            discreteResidualPower =
                returnReduce(discreteResidualPower, sumOp<scalar>());

            directAdvectionFinalPower =
                returnReduce(directAdvectionFinalPower, sumOp<scalar>());

            deferredLagPower =
                returnReduce(deferredLagPower, sumOp<scalar>());

            maxDiscreteResidualPerMass =
                returnReduce
                (
                    maxDiscreteResidualPerMass,
                    maxOp<scalar>()
                );

            maxDeferredLagPerMass =
                returnReduce
                (
                    maxDeferredLagPerMass,
                    maxOp<scalar>()
                );

            const scalar discreteScale =
                max
                (
                    max
                    (
                        mag(sensibleStoragePower),
                        mag(latentStoragePower)
                    ),
                    max
                    (
                        max
                        (
                            mag(effectiveAdvectionPower),
                            mag(conductionPower)
                        ),
                        SMALL
                    )
                );

            const scalar relativeDiscreteResidual =
                discreteResidualPower/discreteScale;

            Info<< "Discrete Eq. (13) energy audit (v11c)" << nl
                << "    sensible storage          = "
                << sensibleStoragePower << " W" << nl
                << "    latent storage            = "
                << latentStoragePower << " W" << nl
                << "    effective advection used  = "
                << effectiveAdvectionPower << " W" << nl
                << "    conduction RHS            = "
                << conductionPower << " W" << nl
                << "    discrete residual         = "
                << discreteResidualPower << " W" << nl
                << "    relative discrete residual= "
                << relativeDiscreteResidual << nl
                << "    max |cell residual|/rho   = "
                << maxDiscreteResidualPerMass << " W/kg" << nl
                << "    direct advection(final T) = "
                << directAdvectionFinalPower << " W" << nl
                << "    deferred-correction lag   = "
                << deferredLagPower << " W" << nl
                << "    max |deferred lag|        = "
                << maxDeferredLagPerMass << " W/kg"
                << endl;
        }

        scalar maxDeltaTIter = 0.0;

        forAll(T_, celli)
        {
            maxDeltaTIter =
                max(maxDeltaTIter, mag(T_[celli] - TBefore[celli]));
        }

        maxDeltaTIter =
            returnReduce(maxDeltaTIter, maxOp<scalar>());

        // First refresh the Lever state with the new temperature so the
        // species equation uses phase fractions/compositions consistent
        // with this thermal correction.
        const scalar maxDeltaFsThermal =
            updatePhaseState(false);

        // Solve mixture solute conservation. This changes Carbon and hence
        // the local liquidus/solidus and phase compositions.
        const scalar maxDeltaCarbon =
            solveSpeciesTransport();

        // Re-close the phase equilibrium after the Carbon correction.
        const bool finalLoop =
            corr == nSolidificationLoops_ - 1;

        const scalar maxDeltaFsSpecies =
            updatePhaseState(finalLoop);

        const scalar maxDeltaFs =
            max(maxDeltaFsThermal, maxDeltaFsSpecies);

        if (finalLoop)
        {
            updateSpeciesDiagnostics(true);
        }

        Info<< "    loop " << corr + 1
            << "/" << nSolidificationLoops_ << nl
            << "        max |delta T|      = "
            << maxDeltaTIter << " K" << nl
            << "        max |CpT adv corr| = "
            << maxEnergyAdvectionCorrection << " W/kg" << nl
            << "        max |delta Carbon| = "
            << maxDeltaCarbon << nl
            << "        max |delta fs|     = "
            << maxDeltaFs
            << endl;
    }
}

void Foam::solvers::continuousCastingMacrosegregation::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::continuousCastingMacrosegregation::momentumTransportCorrector()
{
    viscosity->correct();
    momentumTransport->correct();
}


void Foam::solvers::continuousCastingMacrosegregation::
thermophysicalTransportCorrector()
{}


void Foam::solvers::continuousCastingMacrosegregation::postSolve()
{
    updateBKCDrag(true);

    // Refresh written species diagnostic fields using the final state of
    // the physical time step. The detailed audit is already printed at the
    // final solidification correction.
    updateSpeciesDiagnostics(false);

    // Retain the v10 field-based global balance for comparison with the
    // new v11 discrete-operator audit printed during the final T correction.
    // This does not affect the solution.
    updateEnergyDiagnostics();
}


// ************************************************************************* //
