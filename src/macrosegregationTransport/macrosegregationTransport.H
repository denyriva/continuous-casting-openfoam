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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

Class
    Foam::functionObjects::macrosegregationTransport

Description
    Gate 4C carbon macrosegregation transport function object.

    The conserved mixture carbon mass fraction Carbon is transported with the
    effective carbon mass flux

        phiCarbon =
        [phi - (1-kp)*(1-alpha)*phiS]
        --------------------------------
              alpha + kp*(1-alpha)

    where

        phiS = rho*(Us . Sf)

    and alpha is the liquid fraction.

    The class is based on the Foundation scalarTransport function-object
    architecture and retains the diffusivity modes:

        none
        constant
        viscosity

    Phase concentrations CarbonL and CarbonS are algebraic diagnostics:

        CarbonL = Carbon/[alpha + kp*(1-alpha)]
        CarbonS = kp*CarbonL

SourceFiles
    macrosegregationTransport.C

\*---------------------------------------------------------------------------*/

#ifndef macrosegregationTransport_functionObject_H
#define macrosegregationTransport_functionObject_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "NamedEnum.H"

namespace Foam
{
namespace functionObjects
{

class macrosegregationTransport
:
    public fvMeshFunctionObject
{
public:

    enum class diffusivityType
    {
        none,
        constant,
        viscosity
    };

    static const NamedEnum<diffusivityType, 3> diffusivityTypeNames_;


private:

    word fieldName_;
    word phiName_;
    word rhoName_;
    word alphaName_;
    word UName_;

    diffusivityType diffusivity_;
    scalar D_;
    scalar alphal_;
    scalar alphat_;

    word schemesField_;
    word solverField_;

    scalar partitionCoefficient_;
    vector pullVelocity_;

    volScalarField Carbon_;


    tmp<volScalarField> D() const;


public:

    TypeName("macrosegregationTransport");

    macrosegregationTransport
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    macrosegregationTransport
    (
        const macrosegregationTransport&
    ) = delete;

    virtual ~macrosegregationTransport() = default;

    virtual bool read(const dictionary& dict);

    virtual wordList fields() const;

    virtual bool execute();

    virtual bool write();

    void operator=(const macrosegregationTransport&) = delete;
};

} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //