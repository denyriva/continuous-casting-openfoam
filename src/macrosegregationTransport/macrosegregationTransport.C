/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "macrosegregationTransport.H"

#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "momentumTransportModel.H"
#include "PstreamReduceOps.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(macrosegregationTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        macrosegregationTransport,
        dictionary
    );
}
}


const Foam::NamedEnum
<
    Foam::functionObjects::macrosegregationTransport::diffusivityType,
    3
>
Foam::functionObjects::macrosegregationTransport::diffusivityTypeNames_
{
    "none",
    "constant",
    "viscosity"
};


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::macrosegregationTransport::D() const
{
    switch (diffusivity_)
    {
        case diffusivityType::constant:
        {
            const word Dname(name() + ":D");

            return volScalarField::New
            (
                Dname,
                mesh_,
                dimensionedScalar
                (
                    Dname,
                    dimKinematicViscosity,
                    D_
                )
            );
        }

        case diffusivityType::viscosity:
        {
            const momentumTransportModel& turbulence =
                mesh_.lookupType<momentumTransportModel>();

            return volScalarField::New
            (
                name() + ":D",
                alphal_*turbulence.nu()
              + alphat_*turbulence.nut()
            );
        }

        case diffusivityType::none:
        {
            return volScalarField::New
            (
                name() + ":D",
                mesh_,
                dimensionedScalar(dimKinematicViscosity, 0)
            );
        }
    }

    FatalErrorInFunction
        << "Unhandled diffusivity type"
        << abort(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::functionObjects::macrosegregationTransport::macrosegregationTransport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "Carbon")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    alphaName_
    (
        dict.lookupOrDefault<word>
        (
            "alpha",
            "solidificationMelting1:alpha1"
        )
    ),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    diffusivity_(diffusivityType::none),
    D_(0),
    alphal_(1),
    alphat_(1),
    schemesField_(fieldName_),
    solverField_(fieldName_),
    partitionCoefficient_(0.25),
    pullVelocity_(vector::zero),
    Carbon_
    (
        IOobject
        (
            fieldName_,
            time_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    )
{
    read(dict);
}


bool Foam::functionObjects::macrosegregationTransport::read
(
    const dictionary& dict
)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    fieldName_ = dict.lookupOrDefault<word>("field", fieldName_);
    phiName_ = dict.lookupOrDefault<word>("phi", phiName_);
    rhoName_ = dict.lookupOrDefault<word>("rho", rhoName_);
    alphaName_ = dict.lookupOrDefault<word>("alpha", alphaName_);
    UName_ = dict.lookupOrDefault<word>("U", UName_);

    schemesField_ =
        dict.lookupOrDefault<word>("schemesField", fieldName_);

    solverField_ =
        dict.lookupOrDefault<word>("solverField", fieldName_);

    if (dict.found("diffusivity"))
    {
        diffusivity_ =
            diffusivityTypeNames_.read(dict.lookup("diffusivity"));
    }
    else
    {
        diffusivity_ = diffusivityType::none;
    }

    switch (diffusivity_)
    {
        case diffusivityType::constant:
        {
            D_ = dict.lookup<scalar>("D");
            break;
        }

        case diffusivityType::viscosity:
        {
            alphal_ = dict.lookupOrDefault<scalar>("alphal", 1);
            alphat_ = dict.lookupOrDefault<scalar>("alphat", 1);
            break;
        }

        case diffusivityType::none:
        {
            D_ = 0;
            break;
        }
    }

    partitionCoefficient_ =
        dict.lookupOrDefault<scalar>
        (
            "partitionCoefficient",
            partitionCoefficient_
        );

    if (dict.found("kp"))
    {
        partitionCoefficient_ = dict.lookup<scalar>("kp");
    }

    pullVelocity_ =
        dict.lookupOrDefault<vector>
        (
            "pullVelocity",
            pullVelocity_
        );

    if (dict.found("Us"))
    {
        pullVelocity_ = dict.lookup<vector>("Us");
    }

    if
    (
        partitionCoefficient_ <= 0
     || partitionCoefficient_ > 1
    )
    {
        FatalIOErrorInFunction(dict)
            << "partitionCoefficient must satisfy 0 < kp <= 1. "
            << "Current value: " << partitionCoefficient_
            << exit(FatalIOError);
    }

    Info<< type() << ": field = " << fieldName_ << nl
        << "    phi                    = " << phiName_ << nl
        << "    rho                    = " << rhoName_ << nl
        << "    alpha                  = " << alphaName_ << nl
        << "    partitionCoefficient   = "
        << partitionCoefficient_ << nl
        << "    pullVelocity          = "
        << pullVelocity_ << nl
        << "    diffusivity            = "
        << diffusivityTypeNames_[diffusivity_] << endl;

    return true;
}


Foam::wordList
Foam::functionObjects::macrosegregationTransport::fields() const
{
    return wordList
    ({
        phiName_,
        rhoName_,
        alphaName_
    });
}


bool Foam::functionObjects::macrosegregationTransport::execute()
{
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    const volScalarField& rho =
        mesh_.lookupObject<volScalarField>(rhoName_);

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    if (phi.dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            << "macrosegregationTransport expects a mass flux phi." << nl
            << "    phi field      = " << phiName_ << nl
            << "    dimensions     = " << phi.dimensions() << nl
            << "    expected       = " << dimMass/dimTime
            << abort(FatalError);
    }


    surfaceScalarField alphaF
    (
        IOobject
        (
            name() + ":alphaF",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::interpolate(alpha)
    );

    surfaceScalarField rhoF
    (
        IOobject
        (
            name() + ":rhoF",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        fvc::interpolate(rho)
    );


    surfaceVectorField UsF
    (
        IOobject
        (
            name() + ":UsF",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("Us", dimVelocity, pullVelocity_)
    );

    surfaceScalarField phiS
    (
        IOobject
        (
            name() + ":phiS",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rhoF*(UsF & mesh_.Sf())
    );


    surfaceScalarField denominatorF
    (
        IOobject
        (
            name() + ":denominatorF",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        alphaF
      + partitionCoefficient_*(scalar(1) - alphaF)
    );


    surfaceScalarField phiCarbon
    (
        IOobject
        (
            "phiCarbon",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        (
            phi
          - (scalar(1) - partitionCoefficient_)
           *(scalar(1) - alphaF)
           *phiS
        )
       /denominatorF
    );


    surfaceScalarField phiCarbonVol
    (
        IOobject
        (
            name() + ":phiCarbonVol",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        phiCarbon/rhoF
    );

    const scalarField sumPhiCarbonVol
    (
        fvc::surfaceSum(mag(phiCarbonVol))().primitiveField()
    );

    const scalar maxCarbonCo =
        0.5
       *gMax
        (
            sumPhiCarbonVol
           /mesh_.V().primitiveField()
        )
       *time_.deltaTValue();

    const scalar meanCarbonCo =
        0.5
       *(
            gSum(sumPhiCarbonVol)
           /gSum(mesh_.V().primitiveField())
        )
       *time_.deltaTValue();


    fvModels& fvModels(fvModels::New(mesh_));
    fvConstraints& fvConstraints(fvConstraints::New(mesh_));

    const word divScheme
    (
        "div(phiCarbon," + schemesField_ + ")"
    );

    fvScalarMatrix CarbonEqn
    (
        fvm::ddt(rho, Carbon_)
      + fvm::div(phiCarbon, Carbon_, divScheme)
     ==
        fvModels.source(rho, Carbon_)
    );

    if (diffusivity_ != diffusivityType::none)
    {
        const volScalarField diffusivity(D());

        CarbonEqn -=
            fvm::laplacian
            (
                rho*diffusivity,
                Carbon_,
                "laplacian(" + diffusivity.name()
              + "," + schemesField_ + ")"
            );
    }

    CarbonEqn.relax();
    fvConstraints.constrain(CarbonEqn);

    CarbonEqn.solve(solverField_);

    fvConstraints.constrain(Carbon_);
    Carbon_.correctBoundaryConditions();

    Info<< type() << ":" << nl
        << "    Carbon min/max          = "
        << gMin(Carbon_) << " " << gMax(Carbon_) << nl
        << "    Carbon Courant mean/max = "
        << meanCarbonCo << " " << maxCarbonCo
        << endl;

    return true;
}


bool Foam::functionObjects::macrosegregationTransport::write()
{
    bool ok = Carbon_.write();

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);


    volScalarField CarbonL
    (
        IOobject
        (
            "CarbonL",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Carbon_
       /(
            alpha
          + partitionCoefficient_*(scalar(1) - alpha)
        )
    );

    volScalarField CarbonS
    (
        IOobject
        (
            "CarbonS",
            time_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        partitionCoefficient_*CarbonL
    );

    ok = CarbonL.write() && ok;
    ok = CarbonS.write() && ok;


    if (mesh_.foundObject<volVectorField>(UName_))
    {
        const volVectorField& U =
            mesh_.lookupObject<volVectorField>(UName_);

        const dimensionedVector Us
        (
            "Us",
            dimVelocity,
            pullVelocity_
        );

        volVectorField UCarbon
        (
            IOobject
            (
                "UCarbon",
                time_.name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (
                U
              - (scalar(1) - partitionCoefficient_)
               *(scalar(1) - alpha)
               *Us
            )
           /(
                alpha
              + partitionCoefficient_*(scalar(1) - alpha)
            )
        );

        ok = UCarbon.write() && ok;
    }

    return ok;
}


// ************************************************************************* //