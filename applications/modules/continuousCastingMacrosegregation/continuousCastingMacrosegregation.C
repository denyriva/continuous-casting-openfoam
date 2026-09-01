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
#include "fvcSnGrad.H"
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


Foam::scalar
Foam::solvers::continuousCastingMacrosegregation::maxDeltaT() const
{
    if (!pseudoSteadyEnabled_)
    {
        return basicFluidSolver::maxDeltaT();
    }

    scalar deltaT =
        min
        (
            fvModels().maxDeltaT(),
            min(pseudoDeltaTTarget_, pseudoMaxDeltaT_)
        );

    if (pseudoMaxCo_ < vGreat && CoNum > small)
    {
        deltaT =
            min
            (
                deltaT,
                pseudoMaxCo_/CoNum*runTime.deltaTValue()
            );
    }

    return deltaT;
}


void
Foam::solvers::continuousCastingMacrosegregation::
validatePseudoSteadyProperties() const
{
    if (!pseudoSteadyEnabled_) return;

    if (pseudoMinDeltaT_ <= SMALL || pseudoMaxDeltaT_ < pseudoMinDeltaT_)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "Require 0 < minDeltaT <= maxDeltaT."
            << exit(FatalIOError);
    }

    if (pseudoMaxCo_ <= SMALL)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "maxCo must be > 0." << exit(FatalIOError);
    }

    if (pseudoGrowthFactor_ <= 1.0)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "growthFactor must be > 1." << exit(FatalIOError);
    }

    if (pseudoShrinkFactor_ <= 0.0 || pseudoShrinkFactor_ >= 1.0)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "Require 0 < shrinkFactor < 1." << exit(FatalIOError);
    }

    if (pseudoGrowthThreshold_ <= 0.0 || pseudoGrowthThreshold_ >= 1.0)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "Require 0 < growthThreshold < 1." << exit(FatalIOError);
    }

    if
    (
        pseudoMaxDeltaTStep_ <= SMALL
     || pseudoMaxDeltaCarbonStep_ <= SMALL
     || pseudoMaxDeltaFsStep_ <= SMALL
     || pseudoMaxDeltaUStep_ <= SMALL
    )
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "All per-step change limits must be > 0."
            << exit(FatalIOError);
    }

    if (pseudoReportInterval_ < 1)
    {
        FatalIOErrorInFunction(pseudoSteadyProperties_)
            << "reportInterval must be >= 1."
            << exit(FatalIOError);
    }
}


void
Foam::solvers::continuousCastingMacrosegregation::
updatePseudoSteadyControl()
{
    if (!pseudoSteadyEnabled_) return;

    scalar dTmax=0, dCmax=0, dfsmax=0, dUmax=0;

    const volScalarField& TOld=T_.oldTime();
    const volScalarField& COld=Carbon_.oldTime();
    const volScalarField& fsOld=fs_.oldTime();
    const volVectorField& UOld=U_.oldTime();

    forAll(T_, celli)
    {
        dTmax=max(dTmax, mag(T_[celli]-TOld[celli]));
        dCmax=max(dCmax, mag(Carbon_[celli]-COld[celli]));
        dfsmax=max(dfsmax, mag(fs_[celli]-fsOld[celli]));
        dUmax=max(dUmax, mag(U_[celli]-UOld[celli]));
    }

    dTmax=returnReduce(dTmax,maxOp<scalar>());
    dCmax=returnReduce(dCmax,maxOp<scalar>());
    dfsmax=returnReduce(dfsmax,maxOp<scalar>());
    dUmax=returnReduce(dUmax,maxOp<scalar>());

    const scalar severity =
        max
        (
            dTmax/pseudoMaxDeltaTStep_,
            max
            (
                dCmax/pseudoMaxDeltaCarbonStep_,
                max(dfsmax/pseudoMaxDeltaFsStep_, dUmax/pseudoMaxDeltaUStep_)
            )
        );

    const scalar dt=runTime.deltaTValue();
    word action("hold");

    if (severity > 1.0)
    {
        pseudoDeltaTTarget_=max(pseudoMinDeltaT_,dt*pseudoShrinkFactor_);
        action="shrink";
    }
    else if (severity < pseudoGrowthThreshold_)
    {
        pseudoDeltaTTarget_=
            min(pseudoMaxDeltaT_,max(pseudoMinDeltaT_,dt*pseudoGrowthFactor_));
        action="grow";
    }
    else
    {
        pseudoDeltaTTarget_=min(pseudoMaxDeltaT_,max(pseudoMinDeltaT_,dt));
    }

    if (action=="shrink" || runTime.timeIndex()%pseudoReportInterval_==0)
    {
        Info<< "Pseudo-steady continuation" << nl
            << "    pseudo time         = " << runTime.value() << " s" << nl
            << "    current deltaT      = " << dt << " s" << nl
            << "    next target deltaT  = " << pseudoDeltaTTarget_ << " s" << nl
            << "    max |delta T|       = " << dTmax << " K" << nl
            << "    max |delta Carbon|  = " << dCmax << nl
            << "    max |delta fs|      = " << dfsmax << nl
            << "    max |delta U|       = " << dUmax << " m/s" << nl
            << "    severity            = " << severity << nl
            << "    action              = " << action << endl;
    }
}


bool Foam::solvers::continuousCastingMacrosegregation::diagnosticEnabled
(
    const word& name
) const
{
    if
    (
        !diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "enabled",
            true
        )
    )
    {
        return false;
    }

    const label interval =
        max
        (
            diagnosticsProperties_.lookupOrDefault<label>
            (
                "interval",
                1
            ),
            label(1)
        );

    if (runTime.timeIndex() % interval != 0)
    {
        return false;
    }

    return diagnosticsProperties_.lookupOrDefault<Switch>
    (
        name,
        true
    );
}


Foam::label
Foam::solvers::continuousCastingMacrosegregation::wallCarbonAuditCell() const
{
    const vector target =
        diagnosticsProperties_.lookupOrDefault<vector>
        (
            "wallCarbonTarget",
            vector(-0.0825, 0.0825, 0.5)
        );

    const vectorField& cellCentres = mesh_.C();

    scalar localMinDistSqr = GREAT;
    label localCell = -1;

    forAll(cellCentres, celli)
    {
        const scalar distSqr = magSqr(cellCentres[celli] - target);

        if (distSqr < localMinDistSqr)
        {
            localMinDistSqr = distSqr;
            localCell = celli;
        }
    }

    const scalar globalMinDistSqr =
        returnReduce(localMinDistSqr, minOp<scalar>());

    const scalar matchTolerance =
        SMALL*max(scalar(1), mag(globalMinDistSqr));

    const bool localMatch =
        localCell >= 0
     && mag(localMinDistSqr - globalMinDistSqr) <= matchTolerance;

    // Resolve an unlikely exact tie by selecting the lowest processor rank.
    const label candidateProc =
        localMatch ? Pstream::myProcNo() : Pstream::nProcs();

    const label ownerProc =
        returnReduce(candidateProc, minOp<label>());

    if (localMatch && Pstream::myProcNo() == ownerProc)
    {
        return localCell;
    }

    return -1;
}


void Foam::solvers::continuousCastingMacrosegregation::
updateWallCarbonHistoryDiagnostics() const
{
    // Unlike diagnosticEnabled(), which defaults an unspecified diagnostic
    // to true, this targeted audit is explicitly opt-in.
    if
    (
        !diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "wallCarbonHistory",
            false
        )
     || !diagnosticEnabled("wallCarbonHistory")
    )
    {
        return;
    }

    const label celli = wallCarbonAuditCell();

    if (celli < 0)
    {
        return;
    }

    const scalar fsCell =
        min(max(fs_[celli], scalar(0)), scalar(1));

    const scalar flCell = 1.0 - fsCell;

    const scalar carbonMacrosegregation =
        100.0*(Carbon_[celli] - Carbon0_)/Carbon0_;

    const vector relativeVelocity =
        U_[celli] - solidVelocity_;

    Pout<< "WALL_CARBON_HISTORY"
        << " time=" << runTime.value()
        << " dt=" << runTime.deltaTValue()
        << " proc=" << Pstream::myProcNo()
        << " cell=" << celli
        << " Ccoord=" << mesh_.C()[celli]
        << " T=" << T_[celli]
        << " Carbon=" << Carbon_[celli]
        << " CarbonL=" << CarbonL_[celli]
        << " CarbonS=" << CarbonS_[celli]
        << " CarbonSInterface=" << CarbonSInterface_[celli]
        << " macroseg=" << carbonMacrosegregation
        << " fs=" << fsCell
        << " fl=" << flCell
        << " U=" << U_[celli]
        << " URel=" << relativeVelocity
        << " magURel=" << mag(relativeVelocity)
        << " Dmix=" << speciesDiffusivity_[celli]
        << " liquidAdvectionFactor=" << liquidAdvectionFactor_[celli]
        << " tf=" << localSolidificationTime_[celli]
        << " beta=" << backDiffusionCoefficient_[celli]
        << endl;
}


void Foam::solvers::continuousCastingMacrosegregation::
updateWallCarbonFaceDiagnostics() const
{
    // Explicit opt-in.  This audit compares the actual linearly interpolated
    // liquid-phase diffusion coefficient used at each face of the target
    // cell with a hypothetical harmonic interpolation.  It does not modify
    // the solved species equation.
    if
    (
        !diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "wallCarbonFaceAudit",
            false
        )
     || !diagnosticEnabled("wallCarbonFaceAudit")
    )
    {
        return;
    }

    const label celli = wallCarbonAuditCell();

    const dimensionedScalar DLDim
    (
        "DLWallFaceAudit",
        dimArea/dimTime,
        DL_
    );

    const tmp<volScalarField> tNut = momentumTransport->nut();

    const volScalarField effectiveLiquidDiffusivity
    (
        IOobject
        (
            "wallFaceAuditDlEffTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        DLDim + tNut()/Sct_
    );

    const volScalarField liquidFraction
    (
        IOobject
        (
            "wallFaceAuditLiquidFractionTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        scalar(1) - fs_
    );

    const volScalarField gammaLiquid
    (
        IOobject
        (
            "wallFaceAuditGammaTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        liquidFraction*effectiveLiquidDiffusivity
    );

    // Every rank must execute these coupled surface operations.  In
    // particular, processor-patch interpolation/snGrad exchanges data.
    const tmp<surfaceScalarField> tGammaLinear =
        fvc::interpolate(gammaLiquid);

    const tmp<surfaceScalarField> tSnGradCl =
        fvc::snGrad(CarbonL_);

    const surfaceScalarField& gammaLinear = tGammaLinear();
    const surfaceScalarField& snGradCl = tSnGradCl();

    // Only the rank owning the selected cell prints its faces.
    if (celli < 0)
    {
        return;
    }

    const labelList& cellFaces = mesh_.cells()[celli];
    const labelUList& owner = mesh_.faceOwner();
    const labelUList& neighbour = mesh_.faceNeighbour();

    const scalar cellVolume = mesh_.V()[celli];

    scalar sumLinear = 0.0;
    scalar sumHarmonic = 0.0;

    Pout<< "WALL_CARBON_FACE_BEGIN"
        << " time=" << runTime.value()
        << " proc=" << Pstream::myProcNo()
        << " cell=" << celli
        << " Ccoord=" << mesh_.C()[celli]
        << " fs=" << fs_[celli]
        << " fl=" << 1.0 - fs_[celli]
        << " Carbon=" << Carbon_[celli]
        << " CarbonL=" << CarbonL_[celli]
        << endl;

    forAll(cellFaces, i)
    {
        const label facei = cellFaces[i];
        const bool targetIsOwner = owner[facei] == celli;
        const scalar orientation = targetIsOwner ? 1.0 : -1.0;

        scalar gammaP = gammaLiquid[celli];
        scalar gammaN = gammaP;
        scalar flN = 1.0 - fs_[celli];
        scalar clN = CarbonL_[celli];
        scalar interpolationWeight = 1.0;
        scalar gammaLinearFace = 0.0;
        scalar snGradClFace = 0.0;
        scalar area = 0.0;
        scalar gammaHarmonic = 0.0;
        word faceKind("physicalBoundary");
        word patchName("none");

        if (facei < mesh_.nInternalFaces())
        {
            faceKind = "internal";
            gammaLinearFace = gammaLinear[facei];
            snGradClFace = snGradCl[facei];
            area = mesh_.magSf()[facei];
            gammaHarmonic = gammaLinearFace;

            const label otherCell =
                targetIsOwner ? neighbour[facei] : owner[facei];

            gammaN = gammaLiquid[otherCell];
            flN = 1.0 - fs_[otherCell];
            clN = CarbonL_[otherCell];

            const scalar denom = gammaP - gammaN;

            if (mag(denom) > VSMALL)
            {
                interpolationWeight =
                    (gammaLinearFace - gammaN)/denom;
                interpolationWeight =
                    min(max(interpolationWeight, scalar(0)), scalar(1));
            }
            else
            {
                interpolationWeight = 0.5;
            }

            if (gammaP <= VSMALL || gammaN <= VSMALL)
            {
                gammaHarmonic = 0.0;
            }
            else
            {
                gammaHarmonic =
                    1.0
                   /(
                        interpolationWeight/gammaP
                      + (1.0 - interpolationWeight)/gammaN
                    );
            }
        }
        else
        {
            // OpenFOAM Foundation v14: fvMesh exposes fvBoundaryMesh via
            // boundary(), but not polyMesh::boundaryMesh() directly here.
            // Locate the boundary patch containing this global face index.
            label patchi = -1;

            forAll(mesh_.boundary(), patchj)
            {
                const fvPatch& patch = mesh_.boundary()[patchj];

                if
                (
                    facei >= patch.start()
                 && facei < patch.start() + patch.size()
                )
                {
                    patchi = patchj;
                    break;
                }
            }

            if (patchi < 0)
            {
                FatalErrorInFunction
                    << "Unable to identify boundary patch for face "
                    << facei << exit(FatalError);
            }

            const label patchFacei =
                facei - mesh_.boundary()[patchi].start();

            patchName = mesh_.boundary()[patchi].name();
            gammaLinearFace = gammaLinear.boundaryField()[patchi][patchFacei];
            snGradClFace = snGradCl.boundaryField()[patchi][patchFacei];
            area = mesh_.boundary()[patchi].magSf()[patchFacei];
            gammaHarmonic = gammaLinearFace;

            if (mesh_.boundary()[patchi].coupled())
            {
                faceKind = "processorCoupled";

                const tmp<scalarField> tGammaNeighbour =
                    gammaLiquid.boundaryField()[patchi].patchNeighbourField();

                const tmp<scalarField> tClNeighbour =
                    CarbonL_.boundaryField()[patchi].patchNeighbourField();

                const tmp<scalarField> tFsNeighbour =
                    fs_.boundaryField()[patchi].patchNeighbourField();

                gammaN = tGammaNeighbour()[patchFacei];
                clN = tClNeighbour()[patchFacei];
                flN = 1.0 - tFsNeighbour()[patchFacei];

                const scalar denom = gammaP - gammaN;

                if (mag(denom) > VSMALL)
                {
                    interpolationWeight =
                        (gammaLinearFace - gammaN)/denom;
                    interpolationWeight =
                        min
                        (
                            max(interpolationWeight, scalar(0)),
                            scalar(1)
                        );
                }
                else
                {
                    interpolationWeight = 0.5;
                }

                if (gammaP <= VSMALL || gammaN <= VSMALL)
                {
                    gammaHarmonic = 0.0;
                }
                else
                {
                    gammaHarmonic =
                        1.0
                       /(
                            interpolationWeight/gammaP
                          + (1.0 - interpolationWeight)/gammaN
                        );
                }
            }
            else
            {
                // On a physical boundary there is no neighbouring cell for
                // a harmonic two-cell comparison.  Keep the actual face
                // coefficient for both values; snGrad(Cl) still reveals any
                // boundary-normal liquid diffusion.
                gammaN =
                    gammaLiquid.boundaryField()[patchi][patchFacei];
                flN =
                    1.0 - fs_.boundaryField()[patchi][patchFacei];
                clN =
                    CarbonL_.boundaryField()[patchi][patchFacei];
                gammaHarmonic = gammaLinearFace;
                interpolationWeight = 1.0;
            }
        }

        // This is the face contribution appearing in +div(Gamma grad Cl),
        // oriented outward from the selected target cell.  It is not the
        // Fick flux sign convention (-Gamma grad Cl).
        const scalar operatorFluxLinear =
            orientation*gammaLinearFace*snGradClFace*area;

        const scalar operatorFluxHarmonic =
            orientation*gammaHarmonic*snGradClFace*area;

        sumLinear += operatorFluxLinear;
        sumHarmonic += operatorFluxHarmonic;

        Pout<< "WALL_CARBON_FACE"
            << " time=" << runTime.value()
            << " proc=" << Pstream::myProcNo()
            << " cell=" << celli
            << " face=" << facei
            << " kind=" << faceKind
            << " patch=" << patchName
            << " area=" << area
            << " flP=" << 1.0 - fs_[celli]
            << " flN=" << flN
            << " ClP=" << CarbonL_[celli]
            << " ClN=" << clN
            << " gammaP=" << gammaP
            << " gammaN=" << gammaN
            << " weightP=" << interpolationWeight
            << " gammaLinear=" << gammaLinearFace
            << " gammaHarmonic=" << gammaHarmonic
            << " snGradCl=" << snGradClFace
            << " opFluxLinear=" << operatorFluxLinear
            << " opFluxHarmonic=" << operatorFluxHarmonic
            << endl;
    }

    Pout<< "WALL_CARBON_FACE_SUM"
        << " time=" << runTime.value()
        << " proc=" << Pstream::myProcNo()
        << " cell=" << celli
        << " fs=" << fs_[celli]
        << " fl=" << 1.0 - fs_[celli]
        << " linearTendency=" << sumLinear/cellVolume
        << " harmonicTendency=" << sumHarmonic/cellVolume
        << " ratio="
        <<
        (
            mag(sumLinear) > VSMALL
          ? sumHarmonic/sumLinear
          : scalar(0)
        )
        << endl;
}


void Foam::solvers::continuousCastingMacrosegregation::
updateWallCarbonFluxDiagnostics() const
{
    // Explicit opt-in for this more expensive operator reconstruction.
    if
    (
        !diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "wallCarbonFluxAudit",
            false
        )
     || !diagnosticEnabled("wallCarbonFluxAudit")
    )
    {
        return;
    }

    const label celli = wallCarbonAuditCell();

    // IMPORTANT: do not return on non-owning ranks here.  The fvc::div and
    // fvc::laplacian operations below exchange processor-boundary data, so
    // every MPI rank must execute them in the same order.  Only the final
    // cell extraction/printing is restricted to the owning rank.

    const scalar dt = runTime.deltaTValue();

    if (dt <= SMALL)
    {
        return;
    }

    const dimensionedScalar DLDim
    (
        "DLWallAudit",
        dimArea/dimTime,
        DL_
    );

    const dimensionedScalar DSDim
    (
        "DSWallAudit",
        dimArea/dimTime,
        DS_
    );

    const tmp<volScalarField> tNut =
        momentumTransport->nut();

    const volScalarField effectiveLiquidDiffusivity
    (
        IOobject
        (
            "wallAuditDlEffTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        DLDim + tNut()/Sct_
    );

    const volScalarField liquidFraction
    (
        IOobject
        (
            "wallAuditLiquidFractionTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        scalar(1) - fs_
    );

    const dimensionedVector solidVelocityDim
    (
        "solidVelocityWallAudit",
        dimLength/dimTime,
        solidVelocity_
    );

    const surfaceScalarField solidPhi
    (
        IOobject
        (
            "wallAuditSolidPhiTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        solidVelocityDim & mesh_.Sf()
    );

    const surfaceScalarField relativePhi
    (
        IOobject
        (
            "wallAuditRelativePhiTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phi_ - solidPhi
    );

    // Same Dong-form relative-composition field used by solveSpeciesTransport.
    const volScalarField relativeLiquidComposition
    (
        IOobject
        (
            "wallAuditRelativeCompositionTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fs_*(CarbonL_ - CarbonS_)
    );

    // Reconstruct the final-state physical species equation directly:
    //
    //   dC/dt + div(U C)
    //     = div(fl DlEff grad(Cl))
    //     + div(fs Ds grad(Cs))
    //     - div((U-us) fs(Cl-Cs)).
    //
    // These are physical final-state terms, not the internal Picard
    // decomposition used by phaseLinearized or mixtureCorrection.
    const scalar rDeltaT = 1.0/dt;

    const volScalarField storage
    (
        IOobject
        (
            "wallAuditStorageTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rDeltaT*(Carbon_ - Carbon_.oldTime())
    );

    const tmp<volScalarField> tBulkAdvection =
        fvc::div(phi_, Carbon_, "div(phi,Carbon)");

    const tmp<volScalarField> tLiquidDiffusion =
        fvc::laplacian
        (
            liquidFraction*effectiveLiquidDiffusivity,
            CarbonL_
        );

    const tmp<volScalarField> tSolidDiffusion =
        fvc::laplacian
        (
            fs_*DSDim,
            CarbonS_
        );

    const tmp<volScalarField> tRelativeAdvection =
        fvc::div
        (
            relativePhi,
            relativeLiquidComposition,
            "div(phi,Carbon)"
        );

    const volScalarField& bulkAdvection = tBulkAdvection();
    const volScalarField& liquidDiffusion = tLiquidDiffusion();
    const volScalarField& solidDiffusion = tSolidDiffusion();
    const volScalarField& relativeAdvection = tRelativeAdvection();

    // The parallel operators above must be evaluated on every rank.
    // Only the rank owning the globally selected audit cell reads/prints it.
    if (celli >= 0)
    {
        // Print each transport contribution with the sign it contributes
        // to dC/dt.
        const scalar bulkTendency = -bulkAdvection[celli];
        const scalar liquidDiffusionTendency = liquidDiffusion[celli];
        const scalar solidDiffusionTendency = solidDiffusion[celli];
        const scalar relativeAdvectionTendency = -relativeAdvection[celli];

        const scalar predictedTendency =
            bulkTendency
          + liquidDiffusionTendency
          + solidDiffusionTendency
          + relativeAdvectionTendency;

        const scalar storageTendency = storage[celli];
        const scalar residual = storageTendency - predictedTendency;

        const scalar fsCell =
            min(max(fs_[celli], scalar(0)), scalar(1));

        Pout<< "WALL_CARBON_FLUX"
            << " time=" << runTime.value()
            << " dt=" << dt
            << " proc=" << Pstream::myProcNo()
            << " cell=" << celli
            << " Ccoord=" << mesh_.C()[celli]
            << " Carbon=" << Carbon_[celli]
            << " CarbonOld=" << Carbon_.oldTime()[celli]
            << " dCarbon=" << Carbon_[celli] - Carbon_.oldTime()[celli]
            << " fs=" << fsCell
            << " fl=" << 1.0 - fsCell
            << " magURel=" << mag(U_[celli] - solidVelocity_)
            << " storage=" << storageTendency
            << " bulkAdv=" << bulkTendency
            << " liquidDiff=" << liquidDiffusionTendency
            << " solidDiff=" << solidDiffusionTendency
            << " relativeAdv=" << relativeAdvectionTendency
            << " predicted=" << predictedTendency
            << " residual=" << residual
            << endl;
    }
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

    if (Prt_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Prt must be > 0. Current value: " << Prt_
            << exit(FatalIOError);
    }

    if (Sct_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Sct must be > 0. Current value: " << Sct_
            << exit(FatalIOError);
    }

    if (lambda2_ <= SMALL)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "lambda2 must be > 0. Current value: " << lambda2_
            << exit(FatalIOError);
    }

    if
    (
        microsegregationModel_ != "lever"
     && microsegregationModel_ != "vollerBeckermann"
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "microsegregationModel must be 'lever' or "
            << "'vollerBeckermann'. Current value: "
            << microsegregationModel_
            << exit(FatalIOError);
    }

    if
    (
        speciesDiffusionForm_ != "mixtureCorrection"
     && speciesDiffusionForm_ != "phaseLinearized"
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "speciesDiffusionForm must be 'mixtureCorrection' or "
            << "'phaseLinearized'. Current value: "
            << speciesDiffusionForm_
            << exit(FatalIOError);
    }

    if (vollerBeckermannAlphaC_ < 0)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "vollerBeckermannAlphaC must be >= 0. Current value: "
            << vollerBeckermannAlphaC_
            << exit(FatalIOError);
    }

    if (microsegregationXi_ <= 0)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "microsegregationXi must be > 0. Current value: "
            << microsegregationXi_
            << exit(FatalIOError);
    }

    if (nSolidificationLoops_ < 1)
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "nSolidificationLoops must be >= 1. Current value: "
            << nSolidificationLoops_
            << exit(FatalIOError);
    }

    if
    (
        solidificationIterationMode_ != "fixed"
     && solidificationIterationMode_ != "converged"
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "solidificationIterationMode must be 'fixed' or 'converged'. "
            << "Current value: " << solidificationIterationMode_
            << exit(FatalIOError);
    }

    if
    (
        minSolidificationIterations_ < 1
     || maxSolidificationIterations_ < minSolidificationIterations_
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Require 1 <= minSolidificationIterations <= "
            << "maxSolidificationIterations." << nl
            << "    minSolidificationIterations = "
            << minSolidificationIterations_ << nl
            << "    maxSolidificationIterations = "
            << maxSolidificationIterations_
            << exit(FatalIOError);
    }

    if
    (
        temperatureCouplingTolerance_ <= SMALL
     || carbonCouplingTolerance_ <= SMALL
     || solidFractionCouplingTolerance_ <= SMALL
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "All nonlinear coupling tolerances must be > 0." << nl
            << "    temperatureCouplingTolerance = "
            << temperatureCouplingTolerance_ << nl
            << "    carbonCouplingTolerance = "
            << carbonCouplingTolerance_ << nl
            << "    solidFractionCouplingTolerance = "
            << solidFractionCouplingTolerance_
            << exit(FatalIOError);
    }

    if
    (
        temperatureNonlinearRelaxation_ <= SMALL
     || temperatureNonlinearRelaxation_ > 1.0
     || speciesNonlinearRelaxation_ <= SMALL
     || speciesNonlinearRelaxation_ > 1.0
     || solidFractionNonlinearRelaxation_ <= SMALL
     || solidFractionNonlinearRelaxation_ > 1.0
     || minimumNonlinearRelaxation_ <= SMALL
     || minimumNonlinearRelaxation_ > 1.0
     || minimumNonlinearRelaxation_ > temperatureNonlinearRelaxation_
     || minimumNonlinearRelaxation_ > speciesNonlinearRelaxation_
     || minimumNonlinearRelaxation_ > solidFractionNonlinearRelaxation_
     || nonlinearRelaxationReductionFactor_ <= SMALL
     || nonlinearRelaxationReductionFactor_ >= 1.0
     || nonlinearResidualStallRatio_ <= SMALL
     || nonlinearResidualStallRatio_ > 1.0
     || nonlinearResidualBadIterations_ < 1
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Nonlinear relaxation controls are invalid." << nl
            << "Require initial/minimum relaxations in (0,1], reduction "
            << "factor in (0,1), stall ratio in (0,1], and bad-iteration "
            << "count >= 1." << nl
            << "    temperatureNonlinearRelaxation = "
            << temperatureNonlinearRelaxation_ << nl
            << "    speciesNonlinearRelaxation = "
            << speciesNonlinearRelaxation_ << nl
            << "    solidFractionNonlinearRelaxation = "
            << solidFractionNonlinearRelaxation_ << nl
            << "    minimumNonlinearRelaxation = "
            << minimumNonlinearRelaxation_ << nl
            << "    nonlinearRelaxationReductionFactor = "
            << nonlinearRelaxationReductionFactor_ << nl
            << "    nonlinearResidualStallRatio = "
            << nonlinearResidualStallRatio_ << nl
            << "    nonlinearResidualBadIterations = "
            << nonlinearResidualBadIterations_
            << exit(FatalIOError);
    }

    if
    (
        maxThermophysicalSubcycles_ < 1
     || thermophysicalSubcycleFactor_ < 2
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "Invalid thermophysical subcycling controls." << nl
            << "Require maxThermophysicalSubcycles >= 1 and "
            << "thermophysicalSubcycleFactor >= 2." << nl
            << "    maxThermophysicalSubcycles = "
            << maxThermophysicalSubcycles_ << nl
            << "    thermophysicalSubcycleFactor = "
            << thermophysicalSubcycleFactor_
            << exit(FatalIOError);
    }
}


void Foam::solvers::continuousCastingMacrosegregation::
updateLocalSolidificationTime()
{
    if (microsegregationModel_ != "vollerBeckermann")
    {
        return;
    }

    const scalar dt = runTime.deltaTValue();

    if (dt <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive time step while updating local solidification "
            << "time: " << dt
            << abort(FatalError);
    }

    // Dong et al. neglect remelting. Accumulate physical time only while
    // the cell is currently mushy. Once fs reaches unity the value freezes.
    forAll(fs_, celli)
    {
        const scalar f = fs_[celli];

        if (f > SMALL && f < 1.0 - SMALL)
        {
            localSolidificationTime_[celli] += dt;
        }
    }

    localSolidificationTime_.correctBoundaryConditions();
}


Foam::scalar
Foam::solvers::continuousCastingMacrosegregation::updatePhaseState
(
    const bool report,
    const scalar phaseRelaxation
)
{
    if (phaseRelaxation <= SMALL || phaseRelaxation > 1.0)
    {
        FatalErrorInFunction
            << "phaseRelaxation must satisfy 0 < omega <= 1. Current value: "
            << phaseRelaxation
            << abort(FatalError);
    }

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

    scalar minCarbonSInterface = GREAT;
    scalar maxCarbonSInterface = -GREAT;

    scalar minSolidificationTime = GREAT;
    scalar maxSolidificationTime = -GREAT;

    scalar minBackDiffusion = GREAT;
    scalar maxBackDiffusion = -GREAT;

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

        const scalar fsOldIter = fs_[celli];

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

        // Under-relax the algebraic phase-state update during the inner
        // macro<->micro fixed-point iteration.  The converged fixed point is
        // unchanged; omega_fs=1 recovers the original closure exactly.
        if (phaseRelaxation < 1.0)
        {
            fsCell =
                fsOldIter
              + phaseRelaxation*(fsCell - fsOldIter);

            fsCell = min(max(fsCell, scalar(0)), scalar(1));

            // During a relaxed Picard update the local Jacobian of the
            // relaxed phase map is scaled by the same factor.
            dfsdTCell *= phaseRelaxation;
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

        scalar CarbonLCell = Ccell;
        scalar CarbonSCell = Ccell;
        scalar CarbonSInterfaceCell = kp_*Ccell;
        scalar betaCell = 1.0;

        if (microsegregationModel_ == "lever")
        {
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

            CarbonLCell = Ccell/compositionDenominator;
            CarbonSCell = kp_*CarbonLCell;
            CarbonSInterfaceCell = CarbonSCell;
            betaCell = 1.0;
        }
        else
        {
            // Dong et al., Metals 2017, 7, 209, Eqs. (17)-(23).
            const scalar tf =
                max(localSolidificationTime_[celli], scalar(0));

            const scalar alpha =
                4.0*DS_*tf/sqr(lambda2_);

            const scalar alphaPlus =
                alpha + vollerBeckermannAlphaC_;

            if (alphaPlus <= SMALL)
            {
                FatalErrorInFunction
                    << "Non-positive Voller-Beckermann alpha+." << nl
                    << "    cell      = " << celli << nl
                    << "    alpha     = " << alpha << nl
                    << "    alphaC    = " << vollerBeckermannAlphaC_ << nl
                    << "    tf        = " << tf << " s"
                    << abort(FatalError);
            }

            // Eq. (18):
            // beta = 2*alphaPlus*(1-exp(-1/alphaPlus))
            //        - exp(-1/(2*alphaPlus))
            betaCell =
                2.0*alphaPlus
               *(1.0 - exp(-1.0/alphaPlus))
              - exp(-1.0/(2.0*alphaPlus));

            betaCell =
                min(max(betaCell, scalar(0)), scalar(1));

            if (kp_ >= 1.0 - SMALL)
            {
                CarbonLCell = Ccell;
                CarbonSInterfaceCell = Ccell;
                CarbonSCell = Ccell;
            }
            else
            {
                const scalar betaDenominator =
                    1.0 - betaCell*kp_;

                if (betaDenominator <= SMALL)
                {
                    FatalErrorInFunction
                        << "Invalid Voller-Beckermann denominator." << nl
                        << "    cell    = " << celli << nl
                        << "    beta    = " << betaCell << nl
                        << "    kp      = " << kp_
                        << abort(FatalError);
                }

                const scalar base =
                    1.0 - betaDenominator*fsCell;

                if (base <= SMALL)
                {
                    FatalErrorInFunction
                        << "Invalid Voller-Beckermann concentration base."
                        << nl
                        << "    cell    = " << celli << nl
                        << "    fs      = " << fsCell << nl
                        << "    beta    = " << betaCell << nl
                        << "    kp      = " << kp_ << nl
                        << "    base    = " << base
                        << abort(FatalError);
                }

                const scalar exponent =
                    (kp_ - 1.0)/betaDenominator;

                // Eq. (17): complete diffusion in the liquid means the
                // phase-average liquid concentration equals C_l^*.
                CarbonLCell =
                    Ccell*pow(base, exponent);

                // Eq. (22): local equilibrium applies at the interface.
                CarbonSInterfaceCell =
                    kp_*CarbonLCell;

                if (fsCell > microsegregationXi_)
                {
                    // Eq. (16) gives the phase-average solid
                    // concentration. Use max(fs,xi) only at the
                    // vanishing-solid limit so mixture closure remains
                    // exact for all resolved mushy/solid cells.
                    CarbonSCell =
                        (
                            Ccell
                          - (1.0 - fsCell)*CarbonLCell
                        )
                       /max(fsCell, microsegregationXi_);
                }
                else
                {
                    // At fs=0 the average solid composition is undefined
                    // and inactive in every mixture term.
                    CarbonSCell = CarbonSInterfaceCell;
                }
            }

            if
            (
                CarbonLCell < -SMALL
             || CarbonLCell > 1.0 + SMALL
             || CarbonSCell < -SMALL
             || CarbonSCell > 1.0 + SMALL
            )
            {
                FatalErrorInFunction
                    << "Invalid Voller-Beckermann phase concentration." << nl
                    << "    cell       = " << celli << nl
                    << "    Carbon     = " << Ccell << nl
                    << "    fs         = " << fsCell << nl
                    << "    CarbonL    = " << CarbonLCell << nl
                    << "    CarbonS    = " << CarbonSCell << nl
                    << "    beta       = " << betaCell << nl
                    << "    tf         = " << tf << " s"
                    << abort(FatalError);
            }
        }

        fs_[celli] = fsCell;
        CarbonL_[celli] = CarbonLCell;
        CarbonS_[celli] = CarbonSCell;
        CarbonSInterface_[celli] = CarbonSInterfaceCell;
        backDiffusionCoefficient_[celli] = betaCell;
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

        minCarbonSInterface =
            min(minCarbonSInterface, CarbonSInterfaceCell);
        maxCarbonSInterface =
            max(maxCarbonSInterface, CarbonSInterfaceCell);

        minSolidificationTime =
            min(minSolidificationTime, localSolidificationTime_[celli]);
        maxSolidificationTime =
            max(maxSolidificationTime, localSolidificationTime_[celli]);

        minBackDiffusion = min(minBackDiffusion, betaCell);
        maxBackDiffusion = max(maxBackDiffusion, betaCell);

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
    CarbonSInterface_.correctBoundaryConditions();
    backDiffusionCoefficient_.correctBoundaryConditions();
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

    minCarbonSInterface =
        returnReduce(minCarbonSInterface, minOp<scalar>());
    maxCarbonSInterface =
        returnReduce(maxCarbonSInterface, maxOp<scalar>());

    minSolidificationTime =
        returnReduce(minSolidificationTime, minOp<scalar>());
    maxSolidificationTime =
        returnReduce(maxSolidificationTime, maxOp<scalar>());

    minBackDiffusion =
        returnReduce(minBackDiffusion, minOp<scalar>());
    maxBackDiffusion =
        returnReduce(maxBackDiffusion, maxOp<scalar>());

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

    if (report && diagnosticEnabled("phaseState"))
    {
        Info<< "Phase-state/microsegregation audit" << nl
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
            << "    microsegregation      = "
            << microsegregationModel_ << nl
            << "    CarbonS* min/max      = "
            << minCarbonSInterface << " " << maxCarbonSInterface << nl
            << "    tf min/max            = "
            << minSolidificationTime << " "
            << maxSolidificationTime << " s" << nl
            << "    beta min/max          = "
            << minBackDiffusion << " " << maxBackDiffusion << nl
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

    if (!diagnosticEnabled("buoyancy"))
    {
        return;
    }

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
    // Moving-solid standard BKC:
    //
    //   K^-1 = K0^-1 fsB^2/(1-fsB)^3
    //   K0   = lambda2^2/180
    //   fsB  = min(fs, 0.99)
    //
    // For continuous casting the dendritic skeleton moves with the prescribed
    // solid velocity us. After division of momentum by constant rho, the
    // Darcy acceleration is
    //
    //   a_D = -(mu_l/rho) K^-1 (U-us)
    //       = -(nu_l/K) (U-us).
    //
    // CC-2 introduces only this moving-solid BKC momentum coupling. The
    // species and energy solid-velocity terms remain unchanged until CC-3
    // and CC-4 respectively.

    const scalar K0 = sqr(lambda2_)/180.0;
    const scalar nuLiquid = muLiquid_/rho_;

    scalar minInvK = GREAT;
    scalar maxInvK = -GREAT;
    scalar minDragCoeff = GREAT;
    scalar maxDragCoeff = -GREAT;
    scalar maxDragAcceleration = 0.0;
    scalar maxRelativeUMushy = 0.0;
    scalar maxRelativeUSolidLike = 0.0;
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

        const vector relativeVelocity =
            U_[celli] - solidVelocity_;

        const vector dragAcceleration =
            -dragCoeff*relativeVelocity;

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
            maxRelativeUMushy =
                max(maxRelativeUMushy, mag(relativeVelocity));
        }

        if (fsCell >= 0.99)
        {
            maxRelativeUSolidLike =
                max(maxRelativeUSolidLike, mag(relativeVelocity));
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

    if (!report || !diagnosticEnabled("permeability"))
    {
        return;
    }

    minInvK = returnReduce(minInvK, minOp<scalar>());
    maxInvK = returnReduce(maxInvK, maxOp<scalar>());
    minDragCoeff = returnReduce(minDragCoeff, minOp<scalar>());
    maxDragCoeff = returnReduce(maxDragCoeff, maxOp<scalar>());
    maxDragAcceleration =
        returnReduce(maxDragAcceleration, maxOp<scalar>());
    maxRelativeUMushy =
        returnReduce(maxRelativeUMushy, maxOp<scalar>());
    maxRelativeUSolidLike =
        returnReduce(maxRelativeUSolidLike, maxOp<scalar>());
    integralDragForce =
        returnReduce(integralDragForce, sumOp<scalar>());
    netDragForce =
        returnReduce(netDragForce, sumOp<vector>());

    Info<< "Moving-solid BKC permeability audit" << nl
        << "    solidVelocity          = "
        << solidVelocity_ << " m/s" << nl
        << "    K0                     = "
        << K0 << " m2" << nl
        << "    fsB cap                = 0.99" << nl
        << "    invK min/max           = "
        << minInvK << " " << maxInvK << " 1/m2" << nl
        << "    dragCoeff min/max      = "
        << minDragCoeff << " " << maxDragCoeff << " 1/s" << nl
        << "    max(|a_D|)             = "
        << maxDragAcceleration << " m/s2" << nl
        << "    max(|U-us|) mushy      = "
        << maxRelativeUMushy << " m/s" << nl
        << "    max(|U-us|) fs>=0.99   = "
        << maxRelativeUSolidLike << " m/s" << nl
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
    scalar minNut = GREAT;
    scalar maxNut = -GREAT;
    scalar minDt = GREAT;
    scalar maxDt = -GREAT;
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

    const tmp<volScalarField> tNut =
        momentumTransport->nut();

    const volScalarField& nut =
        tNut();

    forAll(Carbon_, celli)
    {
        const scalar fsCell =
            min(max(fs_[celli], scalar(0)), scalar(1));

        const scalar flCell =
            1.0 - fsCell;

        const scalar nutCell =
            max(nut[celli], scalar(0));

        const scalar DtCell =
            nutCell/Sct_;

        // CC-8: turbulent species diffusion is a liquid-phase contribution.
        // The implicit mixture coefficient is paired with explicit
        // phase-difference corrections in solveSpeciesTransport(), so the
        // converged flux is
        //
        //   fl*(DL + nut/Sct)*grad(Cl) + fs*DS*grad(Cs).
        const scalar DCell =
            fsCell*DS_ + flCell*(DL_ + DtCell);

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

        minNut = min(minNut, nutCell);
        maxNut = max(maxNut, nutCell);

        minDt = min(minDt, DtCell);
        maxDt = max(maxDt, DtCell);

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

    if (!report || !diagnosticEnabled("transport"))
    {
        return;
    }

    minD = returnReduce(minD, minOp<scalar>());
    maxD = returnReduce(maxD, maxOp<scalar>());

    minNut = returnReduce(minNut, minOp<scalar>());
    maxNut = returnReduce(maxNut, maxOp<scalar>());

    minDt = returnReduce(minDt, minOp<scalar>());
    maxDt = returnReduce(maxDt, maxOp<scalar>());

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
        << "    DmixEff min/max         = "
        << minD << " " << maxD << " m2/s" << nl
        << "    nut min/max             = "
        << minNut << " " << maxNut << " m2/s" << nl
        << "    Dt=nut/Sct min/max      = "
        << minDt << " " << maxDt << " m2/s" << nl
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
Foam::solvers::continuousCastingMacrosegregation::solveSpeciesTransport
(
    const bool report,
    const scalar nonlinearRelaxation,
    const scalar subDeltaT,
    const volScalarField& CarbonSubOld
)
{
    // Bennon-Incropera mixture species equation used by the paper, Eq. (19):
    //
    // d(C)/dt + div(U C)
    //   = div(DmixEff grad(C))
    //   + div(fl DlEff grad(Cl-C))
    //   + div(fs Ds grad(Cs-C))
    //   - div((U-us)(Cl-C)),
    //
    // where, for CC-8,
    //
    //   DlEff   = DL + nut/Sct
    //   DmixEff = fs*DS + fl*DlEff.
    //
    // Hence the converged diffusive flux remains phase-separated:
    //
    //   fl*DlEff*grad(Cl) + fs*DS*grad(Cs).
    //
    // For moving solid, us is the prescribed uniform solidVelocity_.
    // The bulk mixture advection remains based on U/phi; only the explicit
    // relative phase-advection term uses U-us.
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
    //   + div(U C)                    : implicit in Carbon
    //   - div((U-us) (Cl-C))           : explicit RHS correction
    //
    // Two diffusion linearisations are selectable at run time:
    //
    // mixtureCorrection (historical)
    // --------------------------------
    //   implicit:  div((fl*DlEff + fs*Ds) grad(C))
    //   explicit:  phase-difference correction terms
    //
    // phaseLinearized (CC10 diagnostic)
    // ----------------------------------
    // Freeze the local phase closure during a Picard iteration:
    //
    //   Cl = ql*C,   Cs = qs*C
    //
    // and split Dong/Bennon-Incropera's direct phase flux
    //
    //   fl*DlEff*grad(Cl) + fs*Ds*grad(Cs)
    //
    // as
    //
    //   (fl*DlEff*ql + fs*Ds*qs)*grad(C)
    //   + fl*DlEff*C*grad(ql)
    //   + fs*Ds*C*grad(qs).
    //
    // The first part is implicit in C; the phase-state-gradient pieces are
    // explicit Picard corrections.  At fixed-point convergence the summed
    // finite-volume diffusion operator recovers the direct phase-separated
    // flux, but without relying on cancellation of Dmix*grad(C) against
    // O(1) phase-difference terms.
    //
    // This switch changes only the numerical linearisation of diffusion,
    // not the continuous species equation.

    updateSpeciesDiagnostics(false);

    if (subDeltaT <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive thermophysical substep: " << subDeltaT
            << abort(FatalError);
    }

    const dimensionedScalar rDeltaTSpecies
    (
        "rDeltaTSpecies",
        dimless/dimTime,
        1.0/subDeltaT
    );


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

    const tmp<volScalarField> tNutSpecies =
        momentumTransport->nut();

    const volScalarField effectiveLiquidDiffusivity
    (
        IOobject
        (
            "effectiveLiquidDiffusivityTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        DLDim + tNutSpecies()/Sct_
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

    // Frozen local phase ratios used by the phaseLinearized Picard form.
    // updatePhaseState() has already supplied CarbonL_/CarbonS_ consistent
    // with the current Carbon_ iterate.
    const dimensionedScalar CarbonLinearisationFloor
    (
        "CarbonLinearisationFloor",
        dimless,
        max(microsegregationXi_, scalar(1e-14))
    );

    const volScalarField carbonForLinearisation
    (
        IOobject
        (
            "carbonForLinearisationTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        max(Carbon_, CarbonLinearisationFloor)
    );

    const volScalarField qLiquid
    (
        IOobject
        (
            "qLiquidSpeciesTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        CarbonL_/carbonForLinearisation
    );

    const volScalarField qSolid
    (
        IOobject
        (
            "qSolidSpeciesTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        CarbonS_/carbonForLinearisation
    );

    const volScalarField phaseLinearizedDiffusivity
    (
        IOobject
        (
            "phaseLinearizedDiffusivityTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        liquidFraction*effectiveLiquidDiffusivity*qLiquid
      + fs_*DSDim*qSolid
    );

    // Uniform solid-velocity face flux and liquid-solid relative flux.
    //
    //   phi_s   = us . Sf
    //   phi_rel = phi - phi_s
    //
    // This changes only the final relative-advection term of Eq. (19).
    // The primary mixture advection term div(U C) continues to use phi_.
    const dimensionedVector solidVelocityDim
    (
        "solidVelocity",
        dimLength/dimTime,
        solidVelocity_
    );

    const surfaceScalarField solidPhi
    (
        IOobject
        (
            "solidVelocityFluxTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        solidVelocityDim & mesh_.Sf()
    );

    const surfaceScalarField relativePhi
    (
        IOobject
        (
            "relativeVelocityFluxTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phi_ - solidPhi
    );

    // Dong-form relative-composition field appearing in the final advective
    // term of Eq. (19).
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
        fs_*(CarbonL_ - CarbonS_)
    );

    const scalarField CarbonBefore(Carbon_.primitiveField());

    // Direct phase-flux Picard corrections.  These are exactly the
    // difference between the frozen direct phase operator and the chosen
    // phaseLinearized implicit Jacobian at the current iterate:
    //
    //   laplacian(fl*DlEff, Cl)
    // - laplacian(fl*DlEff*ql, C)
    //
    // and analogously for the solid phase.
    const tmp<volScalarField> tLiquidPhaseLinearizedCorrection =
        fvc::laplacian
        (
            liquidFraction*effectiveLiquidDiffusivity,
            CarbonL_
        )
      - fvc::laplacian
        (
            liquidFraction*effectiveLiquidDiffusivity*qLiquid,
            Carbon_
        );

    const tmp<volScalarField> tSolidPhaseLinearizedCorrection =
        fvc::laplacian
        (
            fs_*DSDim,
            CarbonS_
        )
      - fvc::laplacian
        (
            fs_*DSDim*qSolid,
            Carbon_
        );

    tmp<fvScalarMatrix> tCarbonEqn;

    if (speciesDiffusionForm_ == "phaseLinearized")
    {
        tCarbonEqn =
        (
            fvm::Sp(rDeltaTSpecies, Carbon_)
          + fvm::div(phi_, Carbon_, "div(phi,Carbon)")
          - fvm::laplacian(phaseLinearizedDiffusivity, Carbon_)
         ==
            rDeltaTSpecies*CarbonSubOld
          + tLiquidPhaseLinearizedCorrection()
          + tSolidPhaseLinearizedCorrection()
          - fvc::div
            (
                relativePhi,
                relativeLiquidComposition,
                "div(phi,Carbon)"
            )
        );
    }
    else
    {
        tCarbonEqn =
        (
            fvm::Sp(rDeltaTSpecies, Carbon_)
          + fvm::div(phi_, Carbon_, "div(phi,Carbon)")
          - fvm::laplacian(speciesDiffusivity_, Carbon_)
         ==
            rDeltaTSpecies*CarbonSubOld
          + fvc::laplacian
            (
                liquidFraction*effectiveLiquidDiffusivity,
                CarbonL_ - Carbon_
            )
          + fvc::laplacian
            (
                fs_*DSDim,
                CarbonS_ - Carbon_
            )
          - fvc::div
            (
                relativePhi,
                relativeLiquidComposition,
                "div(phi,Carbon)"
            )
        );
    }

    fvScalarMatrix& CarbonEqn = tCarbonEqn.ref();

    // CC10c nonlinear equation under-relaxation.
    //
    // This acts on the current Carbon Picard equation, not on any physical
    // diffusion/advection coefficient. omega_C=1 reproduces CC10b exactly.
    const scalar omegaC =
        nonlinearRelaxation > SMALL
      ? nonlinearRelaxation
      : speciesNonlinearRelaxation_;

    CarbonEqn.relax(omegaC);
    CarbonEqn.solve();

    // -----------------------------------------------------------------
    // Pre-phase-state Carbon minimum diagnostic
    // Evaluated immediately after CarbonEqn.solve(), before updatePhaseState().
    // Diagnostic only.
    // -----------------------------------------------------------------

    scalar localMinCarbon = GREAT;
    label localMinCarbonCell = -1;

    forAll(Carbon_, celli)
    {
        if (Carbon_[celli] < localMinCarbon)
        {
            localMinCarbon = Carbon_[celli];
            localMinCarbonCell = celli;
        }
    }

    const scalar globalMinCarbon =
        returnReduce(localMinCarbon, minOp<scalar>());

    // Only evaluate the detailed operator decomposition when Carbon is
    // becoming suspicious, to avoid unnecessary work every Picard solve.
    if (globalMinCarbon < 0.0015)
    {
        const scalar rDeltaTFailure =
            1.0/subDeltaT;

        const volScalarField failureStorage
        (
            IOobject
            (
                "failureStorageTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rDeltaTFailure*(Carbon_ - CarbonSubOld)
        );

        const tmp<volScalarField> tFailureBulkAdvection =
            fvc::div(phi_, Carbon_, "div(phi,Carbon)");

        tmp<volScalarField> tFailureMixtureDiffusion;
        tmp<volScalarField> tFailureLiquidCorrection;
        tmp<volScalarField> tFailureSolidCorrection;

        if (speciesDiffusionForm_ == "phaseLinearized")
        {
            tFailureMixtureDiffusion =
                fvc::laplacian
                (
                    phaseLinearizedDiffusivity,
                    Carbon_
                );

            tFailureLiquidCorrection =
                fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity,
                    CarbonL_
                )
              - fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity*qLiquid,
                    Carbon_
                );

            tFailureSolidCorrection =
                fvc::laplacian
                (
                    fs_*DSDim,
                    CarbonS_
                )
              - fvc::laplacian
                (
                    fs_*DSDim*qSolid,
                    Carbon_
                );
        }
        else
        {
            tFailureMixtureDiffusion =
                fvc::laplacian
                (
                    speciesDiffusivity_,
                    Carbon_
                );

            tFailureLiquidCorrection =
                fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity,
                    CarbonL_ - Carbon_
                );

            tFailureSolidCorrection =
                fvc::laplacian
                (
                    fs_*DSDim,
                    CarbonS_ - Carbon_
                );
        }

        const tmp<volScalarField> tFailureRelativeAdvection =
            fvc::div
            (
                relativePhi,
                relativeLiquidComposition,
                "div(phi,Carbon)"
            );

        const volScalarField& failureBulkAdvection =
            tFailureBulkAdvection();

        const volScalarField& failureMixtureDiffusion =
            tFailureMixtureDiffusion();

        const volScalarField& failureLiquidCorrection =
            tFailureLiquidCorrection();

        const volScalarField& failureSolidCorrection =
            tFailureSolidCorrection();

        const volScalarField& failureRelativeAdvection =
            tFailureRelativeAdvection();

        if
        (
            localMinCarbonCell >= 0
         && mag(localMinCarbon - globalMinCarbon)
            <= SMALL*max(scalar(1), mag(globalMinCarbon))
        )
        {
            const label celli = localMinCarbonCell;

            const scalar CarbonBeforeSolve = CarbonBefore[celli];
            const scalar CarbonAfterSolve = Carbon_[celli];
            const scalar deltaCarbonSolve =
                CarbonAfterSolve - CarbonBeforeSolve;

            scalar liquidCorrectionUsed =
                failureLiquidCorrection[celli];

            scalar solidCorrectionUsed =
                failureSolidCorrection[celli];

            if (speciesDiffusionForm_ == "phaseLinearized")
            {
                liquidCorrectionUsed =
                    tLiquidPhaseLinearizedCorrection()[celli];

                solidCorrectionUsed =
                    tSolidPhaseLinearizedCorrection()[celli];
            }

            const scalar netDiffusionUsed =
                failureMixtureDiffusion[celli]
              + liquidCorrectionUsed
              + solidCorrectionUsed;

            const scalar residual =
                failureStorage[celli]
              + failureBulkAdvection[celli]
              - netDiffusionUsed
              + failureRelativeAdvection[celli];

            Pout<< "CARBON_FAILURE_DIAG"
                << " time=" << runTime.value()
                << " dt=" << subDeltaT
                << " proc=" << Pstream::myProcNo()
                << " cell=" << celli
                << " Ccoord=" << mesh_.C()[celli]
                << " Carbon=" << Carbon_[celli]
                << " CarbonBeforeSolve=" << CarbonBeforeSolve
                << " CarbonAfterSolve=" << CarbonAfterSolve
                << " deltaCarbonSolve=" << deltaCarbonSolve
                << " fs=" << fs_[celli]
                << " CarbonL=" << CarbonL_[celli]
                << " CarbonS=" << CarbonS_[celli]
                << " T=" << T_[celli]
                << " U=" << U_[celli]
                << " storage=" << failureStorage[celli]
                << " bulkAdv=" << failureBulkAdvection[celli]
                << " mixtureDiff=" << failureMixtureDiffusion[celli]
                << " liquidCorrRecomputed=" << failureLiquidCorrection[celli]
                << " liquidCorrUsed=" << liquidCorrectionUsed
                << " solidCorrRecomputed=" << failureSolidCorrection[celli]
                << " solidCorrUsed=" << solidCorrectionUsed
                << " netDiffUsed=" << netDiffusionUsed
                << " relativeAdv=" << failureRelativeAdvection[celli]
                << " residual=" << residual
                << endl;
        }
    }


    // -----------------------------------------------------------------
    // CC-9 open-domain species-conservation audit -- DIAGNOSTIC ONLY
    //
    // This is evaluated immediately after the Carbon solve and before
    // updatePhaseState() changes fs, CarbonL or CarbonS. Therefore the
    // explicit phase-difference diffusion and relative-advection fields are
    // the same frozen fields that entered CarbonEqn.
    //
    // Bringing the implemented equation to the left-hand side gives:
    //
    //   dC/dt
    // + div(phi C)
    // - [active implicit diffusion]
    // - [liquid explicit diffusion correction]
    // - [solid explicit diffusion correction]
    // + div(phiRel,Cl-C)
    // = 0.
    //
    // For mixtureCorrection these are the historical Dmix and phase-
    // difference terms.  For phaseLinearized they are the frozen Picard
    // decomposition of the direct phase-separated flux.
    //
    // Multiplication by rho*V and global summation gives kg/s of solute.
    // Unlike the historical "relative inventory err", this balance remains
    // meaningful in an open caster because inlet/outlet transport is included
    // through the divergence/laplacian operators.
    //
    // The implicit advection/diffusion operators are reconstructed explicitly
    // with the converged Carbon field. With bounded higher-order advection,
    // a small reconstruction lag can remain because the matrix correction was
    // assembled from the pre-solve iterate. The reported balance therefore
    // checks conservation of the solved open-domain equation, while the
    // historical inventory-from-initial-state metric is retained only as a
    // state diagnostic.

    if (report && (runTime.writeTime() || diagnosticEnabled("species")))
    {
        // CC-9 validation cases use first-order Euler time integration.
        // Foundation v14 does not provide fvc::ddt(volScalarField), so
        // reconstruct the Euler storage term explicitly from the current
        // time-step size.
        const scalar rDeltaTSpecies =
            1.0/subDeltaT;

        const volScalarField speciesStorage
        (
            IOobject
            (
                "speciesStorageAuditTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            rDeltaTSpecies*(Carbon_ - CarbonSubOld)
        );

        const tmp<volScalarField> tBulkSpeciesAdvection =
            fvc::div(phi_, Carbon_, "div(phi,Carbon)");

        tmp<volScalarField> tMixtureSpeciesDiffusion;
        tmp<volScalarField> tLiquidSpeciesCorrection;
        tmp<volScalarField> tSolidSpeciesCorrection;

        if (speciesDiffusionForm_ == "phaseLinearized")
        {
            tMixtureSpeciesDiffusion =
                fvc::laplacian(phaseLinearizedDiffusivity, Carbon_);

            tLiquidSpeciesCorrection =
                fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity,
                    CarbonL_
                )
              - fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity*qLiquid,
                    Carbon_
                );

            tSolidSpeciesCorrection =
                fvc::laplacian
                (
                    fs_*DSDim,
                    CarbonS_
                )
              - fvc::laplacian
                (
                    fs_*DSDim*qSolid,
                    Carbon_
                );
        }
        else
        {
            tMixtureSpeciesDiffusion =
                fvc::laplacian(speciesDiffusivity_, Carbon_);

            tLiquidSpeciesCorrection =
                fvc::laplacian
                (
                    liquidFraction*effectiveLiquidDiffusivity,
                    CarbonL_ - Carbon_
                );

            tSolidSpeciesCorrection =
                fvc::laplacian
                (
                    fs_*DSDim,
                    CarbonS_ - Carbon_
                );
        }

        const tmp<volScalarField> tRelativeSpeciesAdvection =
            fvc::div
            (
                relativePhi,
                relativeLiquidComposition,
                "div(phi,Carbon)"
            );

        const volScalarField& bulkSpeciesAdvection =
            tBulkSpeciesAdvection();

        const volScalarField& mixtureSpeciesDiffusion =
            tMixtureSpeciesDiffusion();

        const volScalarField& liquidSpeciesCorrection =
            tLiquidSpeciesCorrection();

        const volScalarField& solidSpeciesCorrection =
            tSolidSpeciesCorrection();

        const volScalarField& relativeSpeciesAdvection =
            tRelativeSpeciesAdvection();

        // -------------------------------------------------------------
        // Spatial species diagnostics -- signed local LHS contributions
        //
        // Units: kg/(m3 s)
        //
        //   storage
        // + bulk advection
        // + mixture diffusion
        // + liquid correction
        // + solid correction
        // + relative advection
        // = residual
        //
        // Diffusive contributions carry a minus sign because the
        // corresponding laplacian operators appear on the RHS of the
        // implemented Carbon equation.
        const dimensionedScalar rhoSpeciesDiagnostic
        (
            "rhoSpeciesDiagnostic",
            dimDensity,
            rho_
        );

        const dimensionedScalar rhoRDeltaTSpeciesDiagnostic
        (
            "rhoRDeltaTSpeciesDiagnostic",
            dimDensity/dimTime,
            rho_*rDeltaTSpecies
        );

        volScalarField carbonStorage
        (
            IOobject
            (
                "carbonStorage",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoRDeltaTSpeciesDiagnostic
           *(Carbon_ - CarbonSubOld)
        );

        volScalarField carbonBulkAdvection
        (
            IOobject
            (
                "carbonBulkAdvection",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoSpeciesDiagnostic*bulkSpeciesAdvection
        );

        // Historical field name retained for A/B plotting.  In
        // phaseLinearized mode this is the active implicit diffusion term,
        // not the historical Dmix operator.
        volScalarField carbonMixtureDiffusion
        (
            IOobject
            (
                "carbonMixtureDiffusion",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rhoSpeciesDiagnostic*mixtureSpeciesDiffusion
        );

        volScalarField carbonLiquidCorrection
        (
            IOobject
            (
                "carbonLiquidCorrection",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rhoSpeciesDiagnostic*liquidSpeciesCorrection
        );

        volScalarField carbonSolidCorrection
        (
            IOobject
            (
                "carbonSolidCorrection",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -rhoSpeciesDiagnostic*solidSpeciesCorrection
        );

        volScalarField carbonRelativeAdvection
        (
            IOobject
            (
                "carbonRelativeAdvection",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoSpeciesDiagnostic*relativeSpeciesAdvection
        );

        volScalarField carbonNetDiffusion
        (
            IOobject
            (
                "carbonNetDiffusion",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            carbonMixtureDiffusion
          + carbonLiquidCorrection
          + carbonSolidCorrection
        );

        volScalarField carbonResidual
        (
            IOobject
            (
                "carbonResidual",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            carbonStorage
          + carbonBulkAdvection
          + carbonNetDiffusion
          + carbonRelativeAdvection
        );

        // Write spatial Carbon diagnostics ONLY at the solver's normal
        // output times.  With writeControl adjustableRunTime, OpenFOAM
        // adjusts deltaT so runTime.writeTime() becomes true at the requested
        // physical output time.  Gating the explicit writes here avoids the
        // historical behaviour of creating processor*/<every time-step>/
        // directories that did not correspond to controlDict output times.
        //
        // The fields are local temporaries, so AUTO_WRITE cannot be relied on;
        // they must still be written explicitly when writeTime() is true.
        if (runTime.writeTime())
        {
            carbonStorage.write();
            carbonBulkAdvection.write();
            carbonMixtureDiffusion.write();
            carbonLiquidCorrection.write();
            carbonSolidCorrection.write();
            carbonRelativeAdvection.write();
            carbonNetDiffusion.write();
            carbonResidual.write();

            if (Pstream::master())
            {
                Info<< "Spatial Carbon diagnostics written at controlDict "
                    << "write time " << runTime.value() << " s" << endl;
            }
        }

        scalar storageRate = 0.0;
        scalar bulkAdvectionRate = 0.0;
        scalar mixtureDiffusionRate = 0.0;
        scalar liquidCorrectionRate = 0.0;
        scalar solidCorrectionRate = 0.0;
        scalar relativeAdvectionRate = 0.0;
        scalar residualRate = 0.0;

        scalar maxResidualPerVolume = 0.0;

        const scalarField& V = mesh_.V();

        forAll(Carbon_, celli)
        {
            const scalar rhoV = rho_*V[celli];

            storageRate +=
                rhoV*speciesStorage[celli];

            bulkAdvectionRate +=
                rhoV*bulkSpeciesAdvection[celli];

            mixtureDiffusionRate +=
                rhoV*mixtureSpeciesDiffusion[celli];

            liquidCorrectionRate +=
                rhoV*liquidSpeciesCorrection[celli];

            solidCorrectionRate +=
                rhoV*solidSpeciesCorrection[celli];

            relativeAdvectionRate +=
                rhoV*relativeSpeciesAdvection[celli];

            const scalar cellResidual =
                speciesStorage[celli]
              + bulkSpeciesAdvection[celli]
              - mixtureSpeciesDiffusion[celli]
              - liquidSpeciesCorrection[celli]
              - solidSpeciesCorrection[celli]
              + relativeSpeciesAdvection[celli];

            residualRate +=
                rhoV*cellResidual;

            maxResidualPerVolume =
                max(maxResidualPerVolume, mag(rho_*cellResidual));
        }

        storageRate =
            returnReduce(storageRate, sumOp<scalar>());

        bulkAdvectionRate =
            returnReduce(bulkAdvectionRate, sumOp<scalar>());

        mixtureDiffusionRate =
            returnReduce(mixtureDiffusionRate, sumOp<scalar>());

        liquidCorrectionRate =
            returnReduce(liquidCorrectionRate, sumOp<scalar>());

        solidCorrectionRate =
            returnReduce(solidCorrectionRate, sumOp<scalar>());

        relativeAdvectionRate =
            returnReduce(relativeAdvectionRate, sumOp<scalar>());

        residualRate =
            returnReduce(residualRate, sumOp<scalar>());

        maxResidualPerVolume =
            returnReduce(maxResidualPerVolume, maxOp<scalar>());

        const scalar balanceScale =
            max
            (
                max
                (
                    mag(storageRate),
                    mag(bulkAdvectionRate)
                ),
                max
                (
                    max
                    (
                        mag(mixtureDiffusionRate),
                        mag(liquidCorrectionRate)
                    ),
                    max
                    (
                        max
                        (
                            mag(solidCorrectionRate),
                            mag(relativeAdvectionRate)
                        ),
                        SMALL
                    )
                )
            );

        const scalar relativeResidual =
            residualRate/balanceScale;

        Info<< "Open-domain species conservation audit (CC-9)" << nl
            << "    storage d(rho*C)/dt      = "
            << storageRate << " kg/s" << nl
            << "    bulk advective term      = "
            << bulkAdvectionRate << " kg/s" << nl
            << "    mixture diffusion term   = "
            << mixtureDiffusionRate << " kg/s" << nl
            << "    liquid correction term   = "
            << liquidCorrectionRate << " kg/s" << nl
            << "    solid correction term    = "
            << solidCorrectionRate << " kg/s" << nl
            << "    relative advective term  = "
            << relativeAdvectionRate << " kg/s" << nl
            << "    signed balance residual  = "
            << residualRate << " kg/s" << nl
            << "    relative balance residual= "
            << relativeResidual << nl
            << "    max |cell residual|      = "
            << maxResidualPerVolume << " kg/(m3 s)"
            << endl;
    }

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


void Foam::solvers::continuousCastingMacrosegregation::updatePhaseCouplingDiagnostics() const
{
    // Diagnostic only: quantify the two chain-rule contributions to the
    // evolution of the Lever-rule solid fraction. No solved equation or
    // field is modified here.
    //
    // In the mushy interval,
    //
    //   fs = (T - Tmelt - m*C)/[(1-kp)(T-Tmelt)]
    //
    // therefore
    //
    //   (dfs/dT)_C = m*C/[(1-kp)(T-Tmelt)^2]
    //   (dfs/dC)_T = -m/[(1-kp)(T-Tmelt)].

    // Keep this diagnostic opt-in. diagnosticEnabled() intentionally
    // defaults unknown switches to true, so require the explicit key too.
    if
    (
        !diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "phaseCoupling",
            false
        )
     || !diagnosticEnabled("phaseCoupling")
    )
    {
        return;
    }

    const scalar dt = runTime.deltaTValue();

    if (dt <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive time step in phase-coupling diagnostic: "
            << dt << abort(FatalError);
    }

    const volScalarField& TOld = T_.oldTime();
    const volScalarField& COld = Carbon_.oldTime();
    const scalarField& V = mesh_.V();

    scalar totalVolume = 0.0;
    scalar mushyVolume = 0.0;
    scalar compositionDominantVolume = 0.0;

    scalar thermalAbsIntegral = 0.0;
    scalar compositionAbsIntegral = 0.0;
    scalar combinedAbsIntegral = 0.0;

    scalar maxThermalRate = 0.0;
    scalar maxCompositionRate = 0.0;
    scalar maxCombinedRate = 0.0;

    forAll(T_, celli)
    {
        totalVolume += V[celli];

        const scalar Tcell = T_[celli];
        const scalar Ccell =
            min(max(Carbon_[celli], scalar(0)), scalar(1));

        const scalar Tliq = Tmelt_ + liquidusSlope_*Ccell;
        const scalar Tsol = Tmelt_ + liquidusSlope_*Ccell/kp_;

        // The analytical Lever derivatives below apply only strictly inside
        // the current mushy interval. Fully liquid/solid cells are skipped.
        if (Tcell >= Tliq || Tcell <= Tsol)
        {
            continue;
        }

        const scalar dTToMelt = Tcell - Tmelt_;
        const scalar oneMinusKp = 1.0 - kp_;

        if (mag(dTToMelt) <= SMALL || mag(oneMinusKp) <= SMALL)
        {
            continue;
        }

        const scalar dfsdT =
            liquidusSlope_*Ccell
           /(oneMinusKp*sqr(dTToMelt));

        const scalar dfsdC =
            -liquidusSlope_
           /(oneMinusKp*dTToMelt);

        const scalar thermalRate =
            dfsdT*(T_[celli] - TOld[celli])/dt;

        const scalar compositionRate =
            dfsdC*(Carbon_[celli] - COld[celli])/dt;

        const scalar combinedRate = thermalRate + compositionRate;

        const scalar aThermal = mag(thermalRate);
        const scalar aComposition = mag(compositionRate);
        const scalar aCombined = mag(combinedRate);

        mushyVolume += V[celli];
        thermalAbsIntegral += aThermal*V[celli];
        compositionAbsIntegral += aComposition*V[celli];
        combinedAbsIntegral += aCombined*V[celli];

        if (aComposition > aThermal)
        {
            compositionDominantVolume += V[celli];
        }

        maxThermalRate = max(maxThermalRate, aThermal);
        maxCompositionRate = max(maxCompositionRate, aComposition);
        maxCombinedRate = max(maxCombinedRate, aCombined);
    }

    totalVolume = returnReduce(totalVolume, sumOp<scalar>());
    mushyVolume = returnReduce(mushyVolume, sumOp<scalar>());
    compositionDominantVolume =
        returnReduce(compositionDominantVolume, sumOp<scalar>());

    thermalAbsIntegral =
        returnReduce(thermalAbsIntegral, sumOp<scalar>());
    compositionAbsIntegral =
        returnReduce(compositionAbsIntegral, sumOp<scalar>());
    combinedAbsIntegral =
        returnReduce(combinedAbsIntegral, sumOp<scalar>());

    maxThermalRate =
        returnReduce(maxThermalRate, maxOp<scalar>());
    maxCompositionRate =
        returnReduce(maxCompositionRate, maxOp<scalar>());
    maxCombinedRate =
        returnReduce(maxCombinedRate, maxOp<scalar>());

    const scalar meanThermalRate =
        mushyVolume > SMALL ? thermalAbsIntegral/mushyVolume : 0.0;

    const scalar meanCompositionRate =
        mushyVolume > SMALL ? compositionAbsIntegral/mushyVolume : 0.0;

    const scalar meanCombinedRate =
        mushyVolume > SMALL ? combinedAbsIntegral/mushyVolume : 0.0;

    const scalar compositionDominantFraction =
        mushyVolume > SMALL
      ? compositionDominantVolume/mushyVolume
      : 0.0;

    const scalar mushyVolumeFraction =
        totalVolume > SMALL ? mushyVolume/totalVolume : 0.0;

    Info<< "Lever-rule T/C phase-coupling diagnostic" << nl
        << "    mushy volume fraction                 = "
        << mushyVolumeFraction << nl
        << "    max |(dfs/dT) dT/dt|                  = "
        << maxThermalRate << " 1/s" << nl
        << "    max |(dfs/dC) dC/dt|                  = "
        << maxCompositionRate << " 1/s" << nl
        << "    max |chain-rule dfs/dt|               = "
        << maxCombinedRate << " 1/s" << nl
        << "    mushy mean |thermal contribution|     = "
        << meanThermalRate << " 1/s" << nl
        << "    mushy mean |composition contribution| = "
        << meanCompositionRate << " 1/s" << nl
        << "    mushy mean |combined contribution|    = "
        << meanCombinedRate << " 1/s" << nl
        << "    mushy volume with |C term|>|T term|   = "
        << compositionDominantFraction
        << endl;
}


void Foam::solvers::continuousCastingMacrosegregation::updateEnergyDiagnostics() const
{
    if (!diagnosticEnabled("energy"))
    {
        return;
    }

    // Global moving-solid Eq. (13) audit.
    //
    // The implemented energy equation is, after division by rho,
    //
    //   d(Cp*T)/dt
    // + div(U Cp T)
    // - d(L fs)/dt
    // - div(L fs U)
    // + div(L fs (U-us))
    // - div((kEff/rho + I_liquid CpL nut/Prt) grad(T))
    // = 0.
    //
    // Hence
    //
    //   E* = integral rho (Cp*T - L*fs) dV,
    //
    // and the open-domain global balance is
    //
    //   dE*/dt
    // + PsensibleAdv
    // - PbulkLatentAdv
    // + PrelativeLatentAdv
    // - Pconduction
    // = 0.
    //
    // The three advective operators are intentionally evaluated separately,
    // matching the solver discretization.  This keeps the audit valid for
    // continuous-casting inlet/outlet boundaries as well as closed cavities.

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

    // ---------------------------------------------------------------
    // Open-boundary advective energy transport
    //
    // Evaluate exactly the same separate finite-volume divergence forms used
    // by the energy equation.  Their global volume integrals are the net
    // outward advective powers through the physical boundary.
    //
    // Sensible:       div(U Cp T)
    // Bulk latent:    div(L fs U)
    // Relative latent div(L fs (U-us))
    //
    // The final residual uses +sensible -bulk +relative, matching TEqn.

    const dimensionedScalar latentHeatDim
    (
        "latentHeat",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),
        latentHeat_
    );

    const dimensionedVector solidVelocityDim
    (
        "solidVelocity",
        dimLength/dimTime,
        solidVelocity_
    );

    const surfaceScalarField solidPhiEnergyAudit
    (
        IOobject
        (
            "solidVelocityFluxEnergyAuditTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        solidVelocityDim & mesh_.Sf()
    );

    const surfaceScalarField relativePhiEnergyAudit
    (
        IOobject
        (
            "relativeVelocityFluxEnergyAuditTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phi_ - solidPhiEnergyAudit
    );

    const volScalarField CpTEnergyAudit
    (
        IOobject
        (
            "CpTEnergyAuditTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        CpMix_*T_
    );

    const tmp<volScalarField> tSensibleAdvectionAudit =
        fvc::div(phi_, CpTEnergyAudit, "div(phi,T)");

    const tmp<volScalarField> tBulkLatentAdvectionAudit =
        latentHeatDim
       *fvc::div
        (
            phi_,
            fs_,
            "div(phi,T)"
        );

    const tmp<volScalarField> tRelativeLatentAdvectionAudit =
        latentHeatDim
       *fvc::div
        (
            relativePhiEnergyAudit,
            fs_,
            "div(phi,T)"
        );

    const volScalarField& sensibleAdvectionAudit =
        tSensibleAdvectionAudit();

    const volScalarField& bulkLatentAdvectionAudit =
        tBulkLatentAdvectionAudit();

    const volScalarField& relativeLatentAdvectionAudit =
        tRelativeLatentAdvectionAudit();

    scalar sensibleAdvectionPower = 0.0;
    scalar bulkLatentAdvectionPower = 0.0;
    scalar relativeLatentAdvectionPower = 0.0;

    forAll(T_, celli)
    {
        sensibleAdvectionPower +=
            rho_*sensibleAdvectionAudit[celli]*V[celli];

        bulkLatentAdvectionPower +=
            rho_*bulkLatentAdvectionAudit[celli]*V[celli];

        relativeLatentAdvectionPower +=
            rho_*relativeLatentAdvectionAudit[celli]*V[celli];
    }

    sensibleAdvectionPower =
        returnReduce(sensibleAdvectionPower, sumOp<scalar>());

    bulkLatentAdvectionPower =
        returnReduce(bulkLatentAdvectionPower, sumOp<scalar>());

    relativeLatentAdvectionPower =
        returnReduce(relativeLatentAdvectionPower, sumOp<scalar>());

    scalar conductiveBoundaryPower = 0.0;

    const tmp<volScalarField> tNutEnergyAudit =
        momentumTransport->nut();

    const volScalarField& nutEnergyAudit =
        tNutEnergyAudit();

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

        const fvPatchScalarField& fsp =
            fs_.boundaryField()[patchi];

        const fvPatchScalarField& nutp =
            nutEnergyAudit.boundaryField()[patchi];

        const scalarField& magSf = patch.magSf();

        const scalarField snGradT(Tp.snGrad());

        forAll(snGradT, facei)
        {
            const scalar fullyLiquidFace =
                (fsp[facei] <= SMALL ? 1.0 : 0.0);

            const scalar kTurbFace =
                rho_*CpLiquid_*fullyLiquidFace
               *max(nutp[facei], scalar(0))/Prt_;

            conductiveBoundaryPower +=
                (kp[facei] + kTurbFace)
               *snGradT[facei]*magSf[facei];
        }
    }

    conductiveBoundaryPower =
        returnReduce(conductiveBoundaryPower, sumOp<scalar>());

    const scalar effectiveEnergyRate =
        (effectiveEnergy - oldEffectiveEnergy)/deltaTValue;

    const scalar balanceResidual =
        effectiveEnergyRate
      + sensibleAdvectionPower
      - bulkLatentAdvectionPower
      + relativeLatentAdvectionPower
      - conductiveBoundaryPower;

    const scalar balanceScale =
        max
        (
            max
            (
                mag(effectiveEnergyRate),
                max
                (
                    mag(sensibleAdvectionPower),
                    max
                    (
                        mag(bulkLatentAdvectionPower),
                        max
                        (
                            mag(relativeLatentAdvectionPower),
                            mag(conductiveBoundaryPower)
                        )
                    )
                )
            ),
            SMALL
        );

    const scalar relativeBalanceResidual =
        balanceResidual/balanceScale;

    Info<< "Post-solve field energy balance" << nl
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
        << "    sensible advective power   = "
        << sensibleAdvectionPower << " W" << nl
        << "    bulk latent advective power= "
        << bulkLatentAdvectionPower << " W" << nl
        << "    relative latent adv. power = "
        << relativeLatentAdvectionPower << " W" << nl
        << "    boundary integral kEffTot*dT/dn = "
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

    diagnosticsProperties_
    (
        IOobject
        (
            "diagnosticsProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),

    pseudoSteadyProperties_
    (
        IOobject
        (
            "pseudoSteadyProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    pseudoSteadyEnabled_(pseudoSteadyProperties_.lookupOrDefault<Switch>("enabled", false)),
    pseudoMinDeltaT_(pseudoSteadyProperties_.lookupOrDefault<scalar>("minDeltaT", 1e-4)),
    pseudoMaxDeltaT_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxDeltaT", 0.1)),
    pseudoMaxCo_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxCo", 2.0)),
    pseudoGrowthFactor_(pseudoSteadyProperties_.lookupOrDefault<scalar>("growthFactor", 1.2)),
    pseudoShrinkFactor_(pseudoSteadyProperties_.lookupOrDefault<scalar>("shrinkFactor", 0.5)),
    pseudoGrowthThreshold_(pseudoSteadyProperties_.lookupOrDefault<scalar>("growthThreshold", 0.35)),
    pseudoMaxDeltaTStep_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxDeltaTStep", 2.0)),
    pseudoMaxDeltaCarbonStep_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxDeltaCarbonStep", 1e-4)),
    pseudoMaxDeltaFsStep_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxDeltaFsStep", 0.02)),
    pseudoMaxDeltaUStep_(pseudoSteadyProperties_.lookupOrDefault<scalar>("maxDeltaUStep", 0.02)),
    pseudoDeltaTTarget_(pseudoSteadyProperties_.lookupOrDefault<scalar>("initialDeltaT", runTime.deltaTValue())),
    pseudoReportInterval_(pseudoSteadyProperties_.lookupOrDefault<label>("reportInterval", 100)),

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
    Prt_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "Prt",
            0.9
        )
    ),
    Sct_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "Sct",
            1.0
        )
    ),
    betaT_(alloyProperties_.lookup<scalar>("betaT")),
    betaC_(alloyProperties_.lookup<scalar>("betaC")),
    TRef_(alloyProperties_.lookup<scalar>("TRef")),
    lambda2_(alloyProperties_.lookup<scalar>("lambda2")),
    microsegregationModel_
    (
        alloyProperties_.lookupOrDefault<word>
        (
            "microsegregationModel",
            "lever"
        )
    ),
    speciesDiffusionForm_
    (
        alloyProperties_.lookupOrDefault<word>
        (
            "speciesDiffusionForm",
            "mixtureCorrection"
        )
    ),
    vollerBeckermannAlphaC_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "vollerBeckermannAlphaC",
            0.1
        )
    ),
    microsegregationXi_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "microsegregationXi",
            1e-12
        )
    ),
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
    solidificationIterationMode_
    (
        alloyProperties_.lookupOrDefault<word>
        (
            "solidificationIterationMode",
            "fixed"
        )
    ),
    minSolidificationIterations_
    (
        alloyProperties_.lookupOrDefault<label>
        (
            "minSolidificationIterations",
            2
        )
    ),
    maxSolidificationIterations_
    (
        alloyProperties_.lookupOrDefault<label>
        (
            "maxSolidificationIterations",
            20
        )
    ),
    temperatureCouplingTolerance_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "temperatureCouplingTolerance",
            1e-3
        )
    ),
    carbonCouplingTolerance_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "carbonCouplingTolerance",
            1e-7
        )
    ),
    solidFractionCouplingTolerance_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "solidFractionCouplingTolerance",
            1e-5
        )
    ),
    temperatureNonlinearRelaxation_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "temperatureNonlinearRelaxation",
            1.0
        )
    ),
    speciesNonlinearRelaxation_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "speciesNonlinearRelaxation",
            1.0
        )
    ),
    solidFractionNonlinearRelaxation_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "solidFractionNonlinearRelaxation",
            1.0
        )
    ),
    adaptiveNonlinearRelaxation_
    (
        alloyProperties_.lookupOrDefault<bool>
        (
            "adaptiveNonlinearRelaxation",
            false
        )
    ),
    minimumNonlinearRelaxation_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "minimumNonlinearRelaxation",
            0.125
        )
    ),
    nonlinearRelaxationReductionFactor_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "nonlinearRelaxationReductionFactor",
            0.5
        )
    ),
    nonlinearResidualStallRatio_
    (
        alloyProperties_.lookupOrDefault<scalar>
        (
            "nonlinearResidualStallRatio",
            0.95
        )
    ),
    nonlinearResidualBadIterations_
    (
        alloyProperties_.lookupOrDefault<label>
        (
            "nonlinearResidualBadIterations",
            2
        )
    ),
    thermophysicalSubcycling_
    (
        alloyProperties_.lookupOrDefault<bool>
        (
            "thermophysicalSubcycling",
            true
        )
    ),
    maxThermophysicalSubcycles_
    (
        alloyProperties_.lookupOrDefault<label>
        (
            "maxThermophysicalSubcycles",
            8
        )
    ),
    thermophysicalSubcycleFactor_
    (
        alloyProperties_.lookupOrDefault<label>
        (
            "thermophysicalSubcycleFactor",
            2
        )
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

    CarbonSInterface_
    (
        IOobject
        (
            "CarbonSInterface",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    localSolidificationTime_
    (
        IOobject
        (
            "localSolidificationTime",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    backDiffusionCoefficient_
    (
        IOobject
        (
            "backDiffusionCoefficient",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0),
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
    CarbonSInterface(CarbonSInterface_),
    localSolidificationTime(localSolidificationTime_),
    backDiffusionCoefficient(backDiffusionCoefficient_),
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
    validatePseudoSteadyProperties();

    if
    (
        pseudoSteadyEnabled_
     && microsegregationModel_ == "vollerBeckermann"
    )
    {
        FatalIOErrorInFunction(alloyProperties_)
            << "vollerBeckermann requires physical-time solidification "
            << "history. Disable pseudoSteady while validating this model."
            << exit(FatalIOError);
    }

    if (pseudoSteadyEnabled_)
    {
        pseudoDeltaTTarget_=
            min(pseudoMaxDeltaT_,max(pseudoMinDeltaT_,pseudoDeltaTTarget_));

        Info<< "continuousCastingMacrosegregation: EXPERIMENTAL pseudo-steady "
            << "continuation ENABLED" << nl
            << "    written runTime is pseudo-time, not physical time" << nl
            << "    initial target deltaT = " << pseudoDeltaTTarget_ << " s" << nl
            << "    min/max deltaT        = " << pseudoMinDeltaT_ << " "
            << pseudoMaxDeltaT_ << " s" << nl
            << "    pseudo maxCo          = " << pseudoMaxCo_ << endl;
    }

    updatePhaseState(true);
    updateSpeciesDiagnostics(true);

    const scalar alphaTLiquid =
        kLiquid_/(rho_*CpLiquid_);

    const scalar alphaTSolid =
        kSolid_/(rho_*CpSolid_);

    Info<< "continuousCastingMacrosegregation: species discretization" << nl
        << "    Eq. (19) advection = term-by-term" << nl
        << "    div(U*C)           = implicit" << nl
        << "    div(U*(Cl-C))      = explicit RHS correction" << nl
        << "    species diffusion  = " << speciesDiffusionForm_ << nl
        << "    phaseLinearized    = direct phase-flux Picard Jacobian"
        << endl;

    Info<< "continuousCastingMacrosegregation: energy discretization" << nl
        << "    Eq. (13) advection = direct CpMix*T product" << nl
        << "    implementation     = implicit split + deferred correction" << nl
        << "    v11c change        = diagnostics only; solved equations unchanged"
        << endl;

    Info<< "continuousCastingMacrosegregation: solid-velocity coupling" << nl
        << "    solidVelocity     = " << solidVelocity_ << " m/s" << nl
        << "    |solidVelocity|   = " << mag(solidVelocity_) << " m/s" << nl
        << "    CC-2 coupling     = moving-solid BKC momentum" << nl
        << "    CC-3 coupling     = moving-solid relative species transport" << nl
        << "    CC-4 coupling     = moving-solid latent-energy transport"
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
        << "    Prt             = " << Prt_ << nl
        << "    Sct             = " << Sct_ << nl
        << "    betaT           = " << betaT_ << " 1/K" << nl
        << "    betaC           = " << betaC_
        << " 1/(solute mass fraction)" << nl
        << "    TRef            = " << TRef_ << " K" << nl
        << "    g               = " << g_.value() << " m/s2" << nl
        << "    lambda2         = " << lambda2_ << " m" << nl
        << "    microsegregationModel = " << microsegregationModel_ << nl
        << "    speciesDiffusionForm  = " << speciesDiffusionForm_ << nl
        << "    VB alphaC       = " << vollerBeckermannAlphaC_ << nl
        << "    solidification iteration mode = "
        << solidificationIterationMode_ << nl
        << "    fixed solidification loops    = "
        << nSolidificationLoops_ << nl
        << "    min/max coupled iterations    = "
        << minSolidificationIterations_ << " "
        << maxSolidificationIterations_ << nl
        << "    coupling tolerances (T,C,fs)  = "
        << temperatureCouplingTolerance_ << " "
        << carbonCouplingTolerance_ << " "
        << solidFractionCouplingTolerance_ << nl
        << "    nonlinear relaxation (T,C,fs) = "
        << temperatureNonlinearRelaxation_ << " "
        << speciesNonlinearRelaxation_ << " "
        << solidFractionNonlinearRelaxation_ << nl
        << "    adaptive nonlinear relaxation = "
        << adaptiveNonlinearRelaxation_ << nl
        << "    adaptive min/reduction/stall   = "
        << minimumNonlinearRelaxation_ << " "
        << nonlinearRelaxationReductionFactor_ << " "
        << nonlinearResidualStallRatio_ << nl
        << "    adaptive bad-iteration count  = "
        << nonlinearResidualBadIterations_ << nl
        << "    thermophysical subcycling     = "
        << thermophysicalSubcycling_ << nl
        << "    max thermophysical subcycles  = "
        << maxThermophysicalSubcycles_ << nl
        << "    thermophysical subcycle factor= "
        << thermophysicalSubcycleFactor_
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
    // Advance the Voller-Beckermann physical-time history from the phase
    // state carried by the previous completed time step.
    updateLocalSolidificationTime();

    if (pseudoSteadyEnabled_)
    {
        U_.oldTime();
        T_.oldTime();
        Carbon_.oldTime();
        fs_.oldTime();
    }

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
    // Moving-solid standard-BKC energy equation following Eq. (13).
    //
    // IMPORTANT CC-4 DISCRETIZATION CHOICE
    // ------------------------------------
    // Keep the two latent convective terms separate:
    //
    //   +div(rho L fs U)
    //   -div(rho L fs (U-us))
    //
    // even though they can be combined algebraically in the continuous
    // equation. This mirrors the validated species treatment: distinct
    // finite-volume advection operators are not collapsed because bounded
    // interpolation can make the discrete operators non-equivalent.
    //
    // After division by constant rho:
    //
    //   d(CpMix*T)/dt + div(U CpMix T)
    //       = div((kEff/rho) grad(T))
    //       + L d(fs)/dt
    //       + div(L fs U)
    //       - div(L fs (U-us)).
    //
    // The local latent-storage term is evaluated directly from the
    // discrete phase-fraction change between physical time levels:
    //
    //   d(fs)/dt ~= (fs^{n+1,k} - fs^n)/dt.
    //
    // This retains composition-driven phase change in the latent balance.
    //
    // For solidVelocity = 0 the two separately discretized latent-advection
    // terms use the same flux and cancel exactly, recovering CC-3.
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

    const scalar physicalDeltaTValue = runTime.deltaTValue();

    if (physicalDeltaTValue <= SMALL)
    {
        FatalErrorInFunction
            << "Non-positive physical time step: " << physicalDeltaTValue
            << abort(FatalError);
    }

    // Preserve previous physical-time algebraic fields while their current
    // values are changed during the inner solidification iterations.
    // These old-time fields are also required by the global Eq. (13)
    // energy-conservation audit.
    CpMix_.oldTime();
    fs_.oldTime();

    // -----------------------------------------------------------------
    // Thermophysical physical-time subcycling checkpoint.
    //
    // The global OpenFOAM Time object is NOT advanced internally.  If the
    // nonlinear T-C-fs solve fails over the full physical deltaT, only the
    // thermophysical fields are restored and the same physical interval is
    // retried as 2, 4, ... equal thermophysical substeps.  Momentum, phi,
    // pressure and turbulence stay frozen exactly as they already do during
    // the nonlinear Picard loop.
    const scalarField TPhysicalStart(T_.primitiveField());
    const scalarField CarbonPhysicalStart(Carbon_.primitiveField());
    const scalarField fsPhysicalStart(fs_.primitiveField());
    const scalarField CarbonLPhysicalStart(CarbonL_.primitiveField());
    const scalarField CarbonSPhysicalStart(CarbonS_.primitiveField());
    const scalarField CarbonSIPhysicalStart(CarbonSInterface_.primitiveField());
    const scalarField betaPhysicalStart(backDiffusionCoefficient_.primitiveField());
    const scalarField speciesDiffPhysicalStart(speciesDiffusivity_.primitiveField());
    const scalarField liquidAdvPhysicalStart(liquidAdvectionFactor_.primitiveField());
    const scalarField macrosegPhysicalStart(macrosegregation_.primitiveField());
    const scalarField CpMixPhysicalStart(CpMix_.primitiveField());
    const scalarField kEffPhysicalStart(kEff_.primitiveField());
    const scalarField dfsdTPhysicalStart(dfsdT_.primitiveField());
    const scalarField CpAppPhysicalStart(CpApp_.primitiveField());

    auto restoreThermophysicalState = [&]()
    {
        T_.primitiveFieldRef() = TPhysicalStart;
        Carbon_.primitiveFieldRef() = CarbonPhysicalStart;
        fs_.primitiveFieldRef() = fsPhysicalStart;
        CarbonL_.primitiveFieldRef() = CarbonLPhysicalStart;
        CarbonS_.primitiveFieldRef() = CarbonSPhysicalStart;
        CarbonSInterface_.primitiveFieldRef() = CarbonSIPhysicalStart;
        backDiffusionCoefficient_.primitiveFieldRef() = betaPhysicalStart;
        speciesDiffusivity_.primitiveFieldRef() = speciesDiffPhysicalStart;
        liquidAdvectionFactor_.primitiveFieldRef() = liquidAdvPhysicalStart;
        macrosegregation_.primitiveFieldRef() = macrosegPhysicalStart;
        CpMix_.primitiveFieldRef() = CpMixPhysicalStart;
        kEff_.primitiveFieldRef() = kEffPhysicalStart;
        dfsdT_.primitiveFieldRef() = dfsdTPhysicalStart;
        CpApp_.primitiveFieldRef() = CpAppPhysicalStart;

        T_.correctBoundaryConditions();
        Carbon_.correctBoundaryConditions();
        fs_.correctBoundaryConditions();
        CarbonL_.correctBoundaryConditions();
        CarbonS_.correctBoundaryConditions();
        CarbonSInterface_.correctBoundaryConditions();
        backDiffusionCoefficient_.correctBoundaryConditions();
        speciesDiffusivity_.correctBoundaryConditions();
        liquidAdvectionFactor_.correctBoundaryConditions();
        macrosegregation_.correctBoundaryConditions();
        CpMix_.correctBoundaryConditions();
        kEff_.correctBoundaryConditions();
        dfsdT_.correctBoundaryConditions();
        CpApp_.correctBoundaryConditions();
    };

    const dimensionedVector solidVelocityDim
    (
        "solidVelocity",
        dimLength/dimTime,
        solidVelocity_
    );

    const surfaceScalarField solidPhi
    (
        IOobject
        (
            "solidVelocityFluxEnergyTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        solidVelocityDim & mesh_.Sf()
    );

    const surfaceScalarField relativePhiEnergy
    (
        IOobject
        (
            "relativeVelocityFluxEnergyTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phi_ - solidPhi
    );

    const tmp<volScalarField> tNutEnergy =
        momentumTransport->nut();

    const volScalarField& nutEnergy =
        tNutEnergy();

    // Dong-faithful turbulent heat switch:
    // turbulent thermal transport is active only in fully liquid cells.
    // In the mushy and solid regions, heat conduction uses phase-weighted
    // molecular conductivity only.
    volScalarField fullyLiquidEnergy
    (
        IOobject
        (
            "fullyLiquidEnergyTmp",
            runTime.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(fullyLiquidEnergy, celli)
    {
        fullyLiquidEnergy[celli] = (fs_[celli] <= SMALL ? 1.0 : 0.0);
    }

    forAll(fullyLiquidEnergy.boundaryField(), patchi)
    {
        scalarField& maskp = fullyLiquidEnergy.boundaryFieldRef()[patchi];
        const fvPatchScalarField& fsp = fs_.boundaryField()[patchi];

        forAll(maskp, facei)
        {
            maskp[facei] = (fsp[facei] <= SMALL ? 1.0 : 0.0);
        }
    }

    const dimensionedScalar CpLiquidDim
    (
        "CpLiquid",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        CpLiquid_
    );

    const bool energyDiagnosticsNow =
        diagnosticEnabled("energy");

    // CC10a nonlinear-coupling diagnostic.
    //
    // Keep this opt-in: diagnosticEnabled() defaults an unspecified named
    // diagnostic to true, whereas CC10a must not change normal log volume
    // unless explicitly requested in diagnosticsProperties.
    const bool couplingDiagnosticsNow =
        diagnosticsProperties_.lookupOrDefault<Switch>
        (
            "coupling",
            false
        )
     && diagnosticEnabled("coupling");

    // CC10b nonlinear-iteration controller.
    //
    // "fixed" reproduces the pre-CC10b nSolidificationLoops behaviour.
    // "converged" permits up to maxSolidificationIterations and stops only
    // when the T, Carbon and fs max-norm changes all satisfy their tolerances.
    // No under-relaxation is introduced in CC10b.
    const bool convergenceControlled =
        solidificationIterationMode_ == "converged";

    const label iterationLimit =
        convergenceControlled
      ? maxSolidificationIterations_
      : nSolidificationLoops_;

    label nThermoSubcycles = 1;

    while (true)
    {
        if (nThermoSubcycles > 1)
        {
            restoreThermophysicalState();
        }

        const scalar subDeltaTValue =
            physicalDeltaTValue/scalar(nThermoSubcycles);

        const dimensionedScalar rDeltaT
        (
            "rDeltaTThermoSubcycle",
            dimless/dimTime,
            1.0/subDeltaTValue
        );

        // IMPORTANT: these are local history carriers, not physical fields.
        // Do NOT clone T_/Carbon_/... directly here.  In particular, cloning
        // T_ also clones its codedMixed patch fields.  codedMixed contains a
        // dynamic-code redirect object tied to the original registered field;
        // cloning it onto an unregistered temporary can invalidate that
        // redirect and segfault when the real T boundary is next updated by
        // fvm::laplacian().  Use plain calculated temporary patch fields and
        // copy only the cell values required by the local transient terms.
        volScalarField TSubOld
        (
            IOobject
            (
                "TSubOldTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", T_.dimensions(), 0.0)
        );
        TSubOld.primitiveFieldRef() = T_.primitiveField();

        volScalarField CarbonSubOld
        (
            IOobject
            (
                "CarbonSubOldTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", Carbon_.dimensions(), 0.0)
        );
        CarbonSubOld.primitiveFieldRef() = Carbon_.primitiveField();

        volScalarField fsSubOld
        (
            IOobject
            (
                "fsSubOldTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", fs_.dimensions(), 0.0)
        );
        fsSubOld.primitiveFieldRef() = fs_.primitiveField();

        volScalarField CpMixSubOld
        (
            IOobject
            (
                "CpMixSubOldTmp",
                runTime.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", CpMix_.dimensions(), 0.0)
        );
        CpMixSubOld.primitiveFieldRef() = CpMix_.primitiveField();

        bool attemptConverged = true;

        Info<< "COUPLING_SUBCYCLE"
            << " time=" << runTime.value()
            << " action=attempt"
            << " nSub=" << nThermoSubcycles
            << " dtPhysical=" << physicalDeltaTValue
            << " dtThermo=" << subDeltaTValue
            << endl;

        for (label thermoSub = 0; thermoSub < nThermoSubcycles; ++thermoSub)
        {
            // Refresh the Dong fully-liquid turbulence switch at the start of
            // every physical thermophysical substep.
            forAll(fullyLiquidEnergy, celli)
            {
                fullyLiquidEnergy[celli] =
                    (fs_[celli] <= SMALL ? 1.0 : 0.0);
            }
            forAll(fullyLiquidEnergy.boundaryField(), patchi)
            {
                scalarField& maskp =
                    fullyLiquidEnergy.boundaryFieldRef()[patchi];
                const fvPatchScalarField& fsp = fs_.boundaryField()[patchi];
                forAll(maskp, facei)
                {
                    maskp[facei] =
                        (fsp[facei] <= SMALL ? 1.0 : 0.0);
                }
            }

            scalar omegaTCurrent = temperatureNonlinearRelaxation_;
            scalar omegaCCurrent = speciesNonlinearRelaxation_;
            scalar omegaFsCurrent = solidFractionNonlinearRelaxation_;
            scalar previousNormalizedResidual = GREAT;
            label consecutiveBadResiduals = 0;
            bool substepConverged = false;

    if (energyDiagnosticsNow)
    {
        Info<< "Solidification/variable-property energy coupling" << nl
        << "    latent storage              = implicit d(fs)/dT + explicit phase residual" << nl
        << "    Eq. (13) sensible advection = direct CpMix*T product" << nl
        << "    CC-4 latent advection       = term-by-term" << nl
        << "    bulk latent flux            = phi" << nl
        << "    relative latent flux        = phi - phiSolid" << nl
        << "    Dong turbulent heat         = fully-liquid only: CpLiquid*nut/Prt" << nl
        << "    Prt                          = " << Prt_ << nl
        << "    nut min/max                  = "
        << gMin(nutEnergy.primitiveField()) << " "
        << gMax(nutEnergy.primitiveField()) << " m2/s" << nl
        << "    alphaT=nut/Prt min/max       = "
        << gMin(nutEnergy.primitiveField())/Prt_ << " "
        << gMax(nutEnergy.primitiveField())/Prt_ << " m2/s" << nl
        << "    solidVelocity               = "
        << solidVelocity_ << " m/s" << nl
        << "    iteration mode              = "
        << solidificationIterationMode_ << nl
        << "    iteration limit             = "
        << iterationLimit
        << endl;
    }

    for (label corr = 0; corr < iterationLimit; ++corr)
    {
        // CC10a: snapshot the coupled iterate before this fixed-point
        // correction. This is diagnostic-only and does not modify the
        // physical fields.
        const scalarField TBeforeIter(T_.primitiveField());
        const scalarField CarbonBeforeIter(Carbon_.primitiveField());
        const scalarField fsBeforeIter(fs_.primitiveField());

        // Hybrid latent-storage linearisation.
        //
        // The total discrete phase change sought over the physical time step
        // is
        //
        //   fs^{n+1} - fs^n.
        //
        // Keep the thermally stiff part implicit through dfs/dT, as in the
        // original apparent-heat-capacity treatment, and place only the
        // phase-change defect not explained by that thermal linearisation on
        // the explicit RHS.  At fixed-point iteration k, define
        //
        //   Rf^k = (fs^k - fs^n)
        //        - (dfs/dT)^k (T^k - T^n).
        //
        // Rf contains the composition-driven contribution and the nonlinear
        // remainder of the Lever-rule closure.  The next thermal correction
        // therefore represents
        //
        //   fs^{k+1} - fs^n
        //     ~= (dfs/dT)^k (T^{k+1} - T^n) + Rf^k.
        //
        // With tLatentCp = -L*dfs/dT (>0 in the mush), the first term is
        // treated implicitly in the matrix and +L*Rf/dt is added to the RHS.
        const tmp<volScalarField> tLatentCp =
            -latentHeatDim*dfsdT_;

        const tmp<volScalarField> tPhaseResidual =
            (fs_ - fsSubOld)
          - dfsdT_*(T_ - TSubOld);

        const tmp<volScalarField> tResidualLatentSource =
            latentHeatDim*rDeltaT*tPhaseResidual();

        // Effective heat-diffusion coefficient after dividing Eq. (13)
        // by rho.
        //
        // Dong Eq. (10): turbulent heat transport is present only in the
        // fully liquid region.  The mush uses the phase-weighted molecular
        // conductivity kEff = fs*kSolid + fl*kLiquid without a nut term.
        //
        //   fully liquid: kLiquid/rho + CpLiquid*nut/Prt
        //   mush:         kEff/rho
        //   solid:        kSolid/rho
        //
        // For laminar models nut=0.
        const tmp<volScalarField> tKbyRho =
            kEff_/rhoDim
          + fullyLiquidEnergy*CpLiquidDim*nutEnergy/Prt_;

        // CC-4 latent-energy advection terms are discretized separately.
        //
        // Bulk term:     div(L fs U)
        // Relative term: div(L fs (U-us))
        //
        // Do not algebraically collapse these operators.
        const tmp<volScalarField> tBulkLatentAdvection =
            latentHeatDim
           *fvc::div
            (
                phi_,
                fs_,
                "div(phi,T)"
            );

        const tmp<volScalarField> tRelativeLatentAdvection =
            latentHeatDim
           *fvc::div
            (
                relativePhiEnergy,
                fs_,
                "div(phi,T)"
            );

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

        if (energyDiagnosticsNow)
        {
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
        }

        fvScalarMatrix TEqn
        (
            fvm::Sp(rDeltaT*CpMix_, T_)
          + fvm::Sp(rDeltaT*tLatentCp(), T_)
          + fvm::div(phiCp, T_, "div(phi,T)")
          - fvm::laplacian(tKbyRho(), T_)
         ==
            rDeltaT*CpMixSubOld*TSubOld
          + rDeltaT*tLatentCp()*TSubOld
          + tResidualLatentSource()
          - tEnergyAdvectionCorrection()
          + tBulkLatentAdvection()
          - tRelativeLatentAdvection()
        );

        // CC10c nonlinear equation under-relaxation.
        //
        // This damps the segregated fixed-point update without changing the
        // converged equation. omega_T=1 reproduces the CC10b matrix exactly.
        TEqn.relax(omegaTCurrent);
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
        // - L*(dfs/dT)*(T-T.oldTime())/dt
        // - L*Rf/dt
        // + div(interpolate(Cp)*phi,T)
        // + [directAdv(old iterate) - splitAdv(old iterate)]
        // - div(L fs U)
        // + div(L fs (U-us))
        // - laplacian(k/rho,T)
        // = 0.
        //
        // Multiplication by rho*V and global summation gives watts.
        //
        // No field used by the solver is modified below.

        const bool finalThermalCorrection =
            corr == iterationLimit - 1;

        if (finalThermalCorrection && energyDiagnosticsNow)
        {
            // Foundation v14 has no fvc::ddt(CpMix_, T_) overload.
            // For this fixed mesh with Euler time integration, the explicit
            // diagnostic equivalent of fvm::Sp(rDeltaT*CpMix_, T_) is
            //
            //   rDeltaT*(Cp*T - Cp.oldTime()*T.oldTime()).
            //
            // This is diagnostic-only and does not alter TEqn.
            const tmp<volScalarField> tDiscreteSensibleStorage =
                rDeltaT
               *(
                    CpMix_*T_
                  - CpMixSubOld*TSubOld
                );

            // Hybrid latent storage exactly matching TEqn:
            //
            //   -L*(dfs/dT)*(T-Told)/dt - L*Rf/dt.
            //
            // Since tLatentCp = -L*dfs/dT, the first contribution is
            // +tLatentCp*(T-Told)/dt on the left-hand side.
            const tmp<volScalarField> tDiscreteLatentStorage =
                rDeltaT*tLatentCp()*(T_ - TSubOld)
              - latentHeatDim*rDeltaT*tPhaseResidual();

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
              - tBulkLatentAdvection()
              + tRelativeLatentAdvection()
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
            scalar bulkLatentAdvectionPower = 0.0;
            scalar relativeLatentAdvectionPower = 0.0;
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

            const volScalarField& bulkLatentAdvection =
                tBulkLatentAdvection();

            const volScalarField& relativeLatentAdvection =
                tRelativeLatentAdvection();

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

                bulkLatentAdvectionPower +=
                    rhoV*bulkLatentAdvection[celli];

                relativeLatentAdvectionPower +=
                    rhoV*relativeLatentAdvection[celli];

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

            bulkLatentAdvectionPower =
                returnReduce
                (
                    bulkLatentAdvectionPower,
                    sumOp<scalar>()
                );

            relativeLatentAdvectionPower =
                returnReduce
                (
                    relativeLatentAdvectionPower,
                    sumOp<scalar>()
                );

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
                            max
                            (
                                max
                                (
                                    mag(bulkLatentAdvectionPower),
                                    mag(relativeLatentAdvectionPower)
                                ),
                                mag(conductionPower)
                            )
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
                << "    bulk latent advection     = "
                << bulkLatentAdvectionPower << " W" << nl
                << "    relative latent advection = "
                << relativeLatentAdvectionPower << " W" << nl
                << "    max |bulk latent advection|     = "
                << gMax(mag(bulkLatentAdvection.primitiveField()))
                << " W/kg" << nl
                << "    max |relative latent advection| = "
                << gMax(mag(relativeLatentAdvection.primitiveField()))
                << " W/kg" << nl
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
        scalar localMaxDeltaTIter = -GREAT;
        label localMaxDeltaTCell = -1;

        forAll(T_, celli)
        {
            const scalar dTCell =
                mag(T_[celli] - TBeforeIter[celli]);

            maxDeltaTIter =
                max(maxDeltaTIter, dTCell);

            if (dTCell > localMaxDeltaTIter)
            {
                localMaxDeltaTIter = dTCell;
                localMaxDeltaTCell = celli;
            }
        }

        maxDeltaTIter =
            returnReduce(maxDeltaTIter, maxOp<scalar>());

        // First refresh the Lever state with the new temperature so the
        // species equation uses phase fractions/compositions consistent
        // with this thermal correction.
        const scalar maxDeltaFsThermal =
            updatePhaseState(false, omegaFsCurrent);

        // Solve mixture solute conservation. This changes Carbon and hence
        // the local liquidus/solidus and phase compositions. On the final
        // correction, audit the open-domain species balance immediately after
        // the solve and before phase-state re-closure changes the frozen fields.
        const bool finalLoop =
            corr == iterationLimit - 1;

        const scalar maxDeltaCarbon =
            solveSpeciesTransport(finalLoop, omegaCCurrent, subDeltaTValue, CarbonSubOld);

        scalar localMaxDeltaCarbon = -GREAT;
        label localMaxDeltaCarbonCell = -1;

        forAll(Carbon_, celli)
        {
            const scalar dCarbonCell =
                mag(Carbon_[celli] - CarbonBeforeIter[celli]);

            if (dCarbonCell > localMaxDeltaCarbon)
            {
                localMaxDeltaCarbon = dCarbonCell;
                localMaxDeltaCarbonCell = celli;
            }
        }

        // Re-close the phase equilibrium after the Carbon correction.
        const scalar maxDeltaFsSpecies =
            updatePhaseState(finalLoop, omegaFsCurrent);

        const scalar maxDeltaFs =
            max(maxDeltaFsThermal, maxDeltaFsSpecies);

        scalar maxDeltaFsIter = 0.0;
        scalar localMaxDeltaFsIter = -GREAT;
        label localMaxDeltaFsCell = -1;

        forAll(fs_, celli)
        {
            const scalar dFsCell =
                mag(fs_[celli] - fsBeforeIter[celli]);

            maxDeltaFsIter =
                max(maxDeltaFsIter, dFsCell);

            if (dFsCell > localMaxDeltaFsIter)
            {
                localMaxDeltaFsIter = dFsCell;
                localMaxDeltaFsCell = celli;
            }
        }

        maxDeltaFsIter =
            returnReduce(maxDeltaFsIter, maxOp<scalar>());

        if (couplingDiagnosticsNow)
        {
            const scalar globalMaxDeltaCarbon =
                returnReduce(localMaxDeltaCarbon, maxOp<scalar>());

            const scalar maxMatchTolT =
                SMALL*max(scalar(1), mag(maxDeltaTIter));

            const scalar maxMatchTolC =
                SMALL*max(scalar(1), mag(globalMaxDeltaCarbon));

            const scalar maxMatchTolFs =
                SMALL*max(scalar(1), mag(maxDeltaFsIter));

            if
            (
                localMaxDeltaTCell >= 0
             && mag(localMaxDeltaTIter - maxDeltaTIter) <= maxMatchTolT
            )
            {
                const label celli = localMaxDeltaTCell;

                Pout<< "COUPLING_MAX_T" << nl
                    << "    time=" << runTime.value() << nl
                    << "    iter=" << corr + 1 << nl
                    << "    proc=" << Pstream::myProcNo() << nl
                    << "    cell=" << celli << nl
                    << "    C=" << mesh_.C()[celli] << nl
                    << "    TBefore=" << TBeforeIter[celli] << nl
                    << "    TAfter=" << T_[celli] << nl
                    << "    CarbonBefore=" << CarbonBeforeIter[celli] << nl
                    << "    CarbonAfter=" << Carbon_[celli] << nl
                    << "    CarbonL=" << CarbonL_[celli] << nl
                    << "    fsBefore=" << fsBeforeIter[celli] << nl
                    << "    fsAfter=" << fs_[celli] << nl
                    << "    dT=" << localMaxDeltaTIter
                    << endl;
            }

            if
            (
                localMaxDeltaCarbonCell >= 0
             && mag(localMaxDeltaCarbon - globalMaxDeltaCarbon)
                <= maxMatchTolC
            )
            {
                const label celli = localMaxDeltaCarbonCell;

                Pout<< "COUPLING_MAX_CARBON" << nl
                    << "    time=" << runTime.value() << nl
                    << "    iter=" << corr + 1 << nl
                    << "    proc=" << Pstream::myProcNo() << nl
                    << "    cell=" << celli << nl
                    << "    C=" << mesh_.C()[celli] << nl
                    << "    TBefore=" << TBeforeIter[celli] << nl
                    << "    TAfter=" << T_[celli] << nl
                    << "    CarbonBefore=" << CarbonBeforeIter[celli] << nl
                    << "    CarbonAfter=" << Carbon_[celli] << nl
                    << "    CarbonL=" << CarbonL_[celli] << nl
                    << "    fsBefore=" << fsBeforeIter[celli] << nl
                    << "    fsAfter=" << fs_[celli] << nl
                    << "    dCarbon=" << localMaxDeltaCarbon
                    << endl;
            }

            if
            (
                localMaxDeltaFsCell >= 0
             && mag(localMaxDeltaFsIter - maxDeltaFsIter) <= maxMatchTolFs
            )
            {
                const label celli = localMaxDeltaFsCell;

                Pout<< "COUPLING_MAX_FS" << nl
                    << "    time=" << runTime.value() << nl
                    << "    iter=" << corr + 1 << nl
                    << "    proc=" << Pstream::myProcNo() << nl
                    << "    cell=" << celli << nl
                    << "    C=" << mesh_.C()[celli] << nl
                    << "    fsBefore=" << fsBeforeIter[celli] << nl
                    << "    fsAfter=" << fs_[celli] << nl
                    << "    T=" << T_[celli] << nl
                    << "    Carbon=" << Carbon_[celli] << nl
                    << "    CarbonL=" << CarbonL_[celli] << nl
                    << "    dfs=" << localMaxDeltaFsIter
                    << endl;
            }
        }

        // Adaptive nonlinear controller.  Normalize each coupled-field
        // change by its requested convergence tolerance and track the worst
        // component.  If the normalized residual fails to improve by at
        // least (1-stallRatio) for several consecutive iterations, reduce
        // all three Picard relaxation factors for the NEXT correction.
        const scalar normalizedResidual =
            max
            (
                max
                (
                    maxDeltaTIter/temperatureCouplingTolerance_,
                    maxDeltaCarbon/carbonCouplingTolerance_
                ),
                maxDeltaFsIter/solidFractionCouplingTolerance_
            );

        if
        (
            adaptiveNonlinearRelaxation_
         && convergenceControlled
         && corr > 0
        )
        {
            if
            (
                normalizedResidual
             >= nonlinearResidualStallRatio_*previousNormalizedResidual
            )
            {
                ++consecutiveBadResiduals;
            }
            else
            {
                consecutiveBadResiduals = 0;
            }

            if
            (
                consecutiveBadResiduals
             >= nonlinearResidualBadIterations_
            )
            {
                const scalar oldOmegaT = omegaTCurrent;
                const scalar oldOmegaC = omegaCCurrent;
                const scalar oldOmegaFs = omegaFsCurrent;

                omegaTCurrent =
                    max
                    (
                        minimumNonlinearRelaxation_,
                        nonlinearRelaxationReductionFactor_*omegaTCurrent
                    );
                omegaCCurrent =
                    max
                    (
                        minimumNonlinearRelaxation_,
                        nonlinearRelaxationReductionFactor_*omegaCCurrent
                    );
                omegaFsCurrent =
                    max
                    (
                        minimumNonlinearRelaxation_,
                        nonlinearRelaxationReductionFactor_*omegaFsCurrent
                    );

                if
                (
                    omegaTCurrent < oldOmegaT - SMALL
                 || omegaCCurrent < oldOmegaC - SMALL
                 || omegaFsCurrent < oldOmegaFs - SMALL
                )
                {
                    Info<< "COUPLING_ADAPT"
                        << " time=" << runTime.value()
                        << " iter=" << corr + 1
                        << " R=" << normalizedResidual
                        << " Rprev=" << previousNormalizedResidual
                        << " omegaT=" << oldOmegaT << "->" << omegaTCurrent
                        << " omegaC=" << oldOmegaC << "->" << omegaCCurrent
                        << " omegaFs=" << oldOmegaFs << "->" << omegaFsCurrent
                        << endl;
                }

                consecutiveBadResiduals = 0;
            }
        }

        previousNormalizedResidual = normalizedResidual;

        const bool couplingConverged =
            convergenceControlled
         && corr + 1 >= minSolidificationIterations_
         && maxDeltaTIter <= temperatureCouplingTolerance_
         && maxDeltaCarbon <= carbonCouplingTolerance_
         && maxDeltaFsIter <= solidFractionCouplingTolerance_;

        if (couplingDiagnosticsNow)
        {
            Info<< "COUPLING_DIAG"
                << " time=" << runTime.value()
                << " iter=" << corr + 1
                << " limit=" << iterationLimit
                << " mode=" << solidificationIterationMode_
                << " omegaT=" << omegaTCurrent
                << " omegaC=" << omegaCCurrent
                << " omegaFs=" << omegaFsCurrent
                << " dT=" << maxDeltaTIter
                << " dCarbon=" << maxDeltaCarbon
                << " dfs=" << maxDeltaFsIter
                << " dfsThermal=" << maxDeltaFsThermal
                << " dfsSpecies=" << maxDeltaFsSpecies
                << " converged=" << couplingConverged
                << endl;
        }

        if (finalLoop)
        {
            updateSpeciesDiagnostics(true);
        }

        if (energyDiagnosticsNow)
        {
            Info<< "    loop " << corr + 1
                << "/" << iterationLimit << nl
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

        if (couplingConverged)
        {
            // The detailed matrix-level energy/species audits are historically
            // tied to a pre-known final correction. Do not re-solve or
            // re-close fields merely for reporting after an early exit.
            // Final-state field diagnostics are refreshed in postSolve().
            substepConverged = true;

            Info<< "COUPLING_CONTROL"
                << " time=" << runTime.value()
                << " status=converged"
                << " thermoSub=" << thermoSub + 1 << "/" << nThermoSubcycles
                << " dtThermo=" << subDeltaTValue
                << " iterations=" << corr + 1
                << " dT=" << maxDeltaTIter
                << " dCarbon=" << maxDeltaCarbon
                << " dfs=" << maxDeltaFsIter
                << endl;

            break;
        }

        if (convergenceControlled && finalLoop)
        {
            Info<< "COUPLING_SUBCYCLE"
                << " time=" << runTime.value()
                << " action=substep-failed"
                << " thermoSub=" << thermoSub + 1 << "/" << nThermoSubcycles
                << " nSub=" << nThermoSubcycles
                << " dtThermo=" << subDeltaTValue
                << " iterations=" << iterationLimit
                << " dT=" << maxDeltaTIter
                << " dCarbon=" << maxDeltaCarbon
                << " dfs=" << maxDeltaFsIter
                << endl;
        }
    }

            // Legacy fixed-loop mode has no convergence test by definition.
            // Preserve its historical behavior: completion of the requested
            // number of corrections constitutes a successful substep.
            if (!convergenceControlled)
            {
                substepConverged = true;
            }

            if (!substepConverged)
            {
                attemptConverged = false;
                break;
            }

            // Advance ONLY the local thermophysical time history after a
            // successfully converged substep. The OpenFOAM global oldTime()
            // fields remain untouched until the physical timestep completes.
            TSubOld = T_;
            CarbonSubOld = Carbon_;
            fsSubOld = fs_;
            CpMixSubOld = CpMix_;
        }

        if (attemptConverged)
        {
            Info<< "COUPLING_SUBCYCLE"
                << " time=" << runTime.value()
                << " action=accepted"
                << " nSub=" << nThermoSubcycles
                << " dtPhysical=" << physicalDeltaTValue
                << " dtThermo=" << subDeltaTValue
                << endl;
            return;
        }

        const label nextSubcycles =
            nThermoSubcycles*thermophysicalSubcycleFactor_;

        if
        (
            !thermophysicalSubcycling_
         || nextSubcycles > maxThermophysicalSubcycles_
        )
        {
            restoreThermophysicalState();

            FatalErrorInFunction
                << "Thermophysical nonlinear coupling failed and the "
                << "physical timestep will NOT be accepted." << nl
                << "    time                  = " << runTime.value() << nl
                << "    physical deltaT       = " << physicalDeltaTValue << nl
                << "    attempted subcycles   = " << nThermoSubcycles << nl
                << "    max allowed subcycles = "
                << maxThermophysicalSubcycles_ << nl
                << "Increase maxThermophysicalSubcycles only after checking "
                << "the coupling diagnostics."
                << abort(FatalError);
        }

        Info<< "COUPLING_SUBCYCLE"
            << " time=" << runTime.value()
            << " action=restore-and-retry"
            << " failedNSub=" << nThermoSubcycles
            << " retryNSub=" << nextSubcycles
            << " retryDtThermo="
            << physicalDeltaTValue/scalar(nextSubcycles)
            << endl;

        nThermoSubcycles = nextSubcycles;
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
    // Per-time-step flow/Courant diagnostic
    const scalar flowDiagDt = runTime.deltaTValue();

    const scalarField flowDiagSumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()
    );


    scalar localMaxCo = -GREAT;
    scalar localMaxU = -GREAT;
    scalar localMaxURel = -GREAT;

    label localMaxCoCell = -1;
    label localMaxUCell = -1;
    label localMaxURelCell = -1;

    forAll(U_, celli)
    {
        const scalar cellCo =
            0.5*flowDiagDt*flowDiagSumPhi[celli]/mesh_.V()[celli];

        const scalar magUCell = mag(U_[celli]);
        const scalar magURelCell = mag(U_[celli] - solidVelocity_);

        if (cellCo > localMaxCo)
        {
            localMaxCo = cellCo;
            localMaxCoCell = celli;
        }

        if (magUCell > localMaxU)
        {
            localMaxU = magUCell;
            localMaxUCell = celli;
        }

        if (magURelCell > localMaxURel)
        {
            localMaxURel = magURelCell;
            localMaxURelCell = celli;
        }
    }

    const scalar globalMaxCo =
        returnReduce(localMaxCo, maxOp<scalar>());

    const scalar globalMaxU =
        returnReduce(localMaxU, maxOp<scalar>());

    const scalar globalMaxURel =
        returnReduce(localMaxURel, maxOp<scalar>());

    if
    (
        localMaxCoCell >= 0
     && mag(localMaxCo - globalMaxCo)
        <= SMALL*max(scalar(1), mag(globalMaxCo))
    )
    {
        Pout<< "FLOW_DIAG CoMax"
            << " time=" << runTime.value()
            << " dt=" << flowDiagDt
            << " Co=" << localMaxCo
            << " proc=" << Pstream::myProcNo()
            << " cell=" << localMaxCoCell
            << " C=" << mesh_.C()[localMaxCoCell]
            << " U=" << U_[localMaxCoCell]
            << " |U|=" << mag(U_[localMaxCoCell])
            << " |U-us|="
            << mag(U_[localMaxCoCell] - solidVelocity_)
            << " fs=" << fs_[localMaxCoCell]
            << " Carbon=" << Carbon_[localMaxCoCell]
            << endl;
    }

    if
    (
        localMaxUCell >= 0
     && mag(localMaxU - globalMaxU)
        <= SMALL*max(scalar(1), mag(globalMaxU))
    )
    {
        Pout<< "FLOW_DIAG UMax"
            << " time=" << runTime.value()
            << " |U|=" << localMaxU
            << " proc=" << Pstream::myProcNo()
            << " cell=" << localMaxUCell
            << " C=" << mesh_.C()[localMaxUCell]
            << " U=" << U_[localMaxUCell]
            << " fs=" << fs_[localMaxUCell]
            << " Carbon=" << Carbon_[localMaxUCell]
            << endl;
    }

    if
    (
        localMaxURelCell >= 0
     && mag(localMaxURel - globalMaxURel)
        <= SMALL*max(scalar(1), mag(globalMaxURel))
    )
    {
        Pout<< "FLOW_DIAG URelMax"
            << " time=" << runTime.value()
            << " |U-us|=" << localMaxURel
            << " proc=" << Pstream::myProcNo()
            << " cell=" << localMaxURelCell
            << " C=" << mesh_.C()[localMaxURelCell]
            << " U=" << U_[localMaxURelCell]
            << " fs=" << fs_[localMaxURelCell]
            << " Carbon=" << Carbon_[localMaxURelCell]
            << endl;
    }


    // Refresh written species diagnostic fields using the final state of
    // the physical time step. The detailed audit is already printed at the
    // final solidification correction.
    updateSpeciesDiagnostics(false);

    // Targeted wall-adjacent Carbon diagnostics. These reconstruct the
    // accepted physical-time-step state only and never modify the solution.
    updateWallCarbonHistoryDiagnostics();
    updateWallCarbonFluxDiagnostics();
    updateWallCarbonFaceDiagnostics();

    // Retain the v10 field-based global balance for comparison with the
    // new v11 discrete-operator audit printed during the final T correction.
    // This does not affect the solution.
    updateEnergyDiagnostics();

    // Diagnostic-only decomposition of Lever-rule phase evolution into
    // temperature-driven and composition-driven contributions.
    updatePhaseCouplingDiagnostics();

    // Adapt only the NEXT pseudo step. There is deliberately no rollback.
    updatePseudoSteadyControl();
}


// ************************************************************************* //
