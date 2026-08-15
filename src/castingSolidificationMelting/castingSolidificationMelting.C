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

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "castingSolidificationMelting.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(castingSolidificationMelting, 0);

    addToRunTimeSelectionTable(fvModel, castingSolidificationMelting, dictionary);
}
}


const Foam::NamedEnum<Foam::fv::castingSolidificationMelting::thermoMode, 2>
Foam::fv::castingSolidificationMelting::thermoModeTypeNames_
{
    "thermo",
    "lookup"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::castingSolidificationMelting::readCoeffs(const dictionary& dict)
{
    Tsol_ = dict.lookup<scalar>("Tsol");
    Tliq_ = dict.lookupOrDefault<scalar>("Tliq", Tsol_);
    alpha1e_ = dict.lookupOrDefault<scalar>("alpha1e", 0.0);

    compositionDependentLiquidus_ =
        dict.lookupOrDefault<bool>("compositionDependentLiquidus", false);

    CarbonName_ = dict.lookupOrDefault<word>("Carbon", "Carbon");

    if (compositionDependentLiquidus_)
    {
        partitionCoefficient_ =
            dict.lookup<scalar>("partitionCoefficient");
        CarbonRef_ = dict.lookup<scalar>("CarbonRef");
        liquidusSlope_ = dict.lookup<scalar>("liquidusSlope");

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

        if (CarbonRef_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "CarbonRef must be non-negative. Current value: "
                << CarbonRef_ << exit(FatalIOError);
        }

        if (mag(liquidusSlope_) < SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "liquidusSlope must be non-zero when "
                << "compositionDependentLiquidus is enabled."
                << exit(FatalIOError);
        }

        if (mag(alpha1e_) > SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Gate 5A composition-dependent closure currently requires "
                << "alpha1e = 0. Current value: " << alpha1e_
                << exit(FatalIOError);
        }
    }
    else
    {
        // Safe inactive defaults. These values are not used by legacy mode.
        partitionCoefficient_ = 1.0;
        CarbonRef_ = 0.0;
        liquidusSlope_ = 0.0;
    }

    L_ = dict.lookup<scalar>("L");

    relax_ = dict.lookupOrDefault("relax", 0.9);

    mode_ = thermoModeTypeNames_.read(dict.lookup("thermoMode"));

    rhoRef_ = dict.lookup<scalar>("rhoRef");
    TName_ = dict.lookupOrDefault<word>("T", "T");
    CpName_ = dict.lookupOrDefault<word>("Cp", "Cp");
    UName_ = dict.lookupOrDefault<word>("U", "U");
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    Cu_ = dict.lookupOrDefault<scalar>("Cu", 100000);
    q_ = dict.lookupOrDefault("q", 0.001);
    pullVelocity_ =
        dict.lookupOrDefault<vector>("pullVelocity", vector::zero);

    Info<< "castingSolidificationMelting: pullVelocity = "
        << pullVelocity_ << endl;

    if (compositionDependentLiquidus_)
    {
        Info<< "castingSolidificationMelting: Gate 5A composition-dependent "
            << "phase equilibrium enabled" << nl
            << "    Carbon field          = " << CarbonName_ << nl
            << "    partitionCoefficient  = " << partitionCoefficient_ << nl
            << "    CarbonRef             = " << CarbonRef_ << nl
            << "    TliqRef               = " << Tliq_ << " K" << nl
            << "    liquidusSlope         = " << liquidusSlope_
            << " K/(mass fraction)" << endl;
    }
    else
    {
        Info<< "castingSolidificationMelting: fixed Tsol/Tliq phase "
            << "equilibrium enabled" << endl;
    }

    beta_ = dict.lookup<scalar>("beta");

    if (mode_ == thermoMode::lookup)
    {
        CpRef_ = dict.lookup<scalar>("CpRef");
    }

    if (!mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        g_ = dict.lookup("g");
    }
}


Foam::tmp<Foam::volScalarField>
Foam::fv::castingSolidificationMelting::Cp() const
{
    switch (mode_)
    {
        case thermoMode::thermo:
        {
            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            return thermo.Cp();
            break;
        }
        case thermoMode::lookup:
        {
            if (CpName_ == "CpRef")
            {
                return volScalarField::New
                (
                    name() + ":Cp",
                    mesh(),
                    dimensionedScalar
                    (
                        dimEnergy/dimMass/dimTemperature,
                        CpRef_
                    ),
                    extrapolatedCalculatedFvPatchScalarField::typeName
                );
            }
            else
            {
                return mesh().lookupObject<volScalarField>(CpName_);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled thermo mode: " << thermoModeTypeNames_[mode_]
                << abort(FatalError);
        }
    }

    return tmp<volScalarField>(nullptr);
}


Foam::vector Foam::fv::castingSolidificationMelting::g() const
{
    if (mesh().foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& value =
            mesh().lookupObject<uniformDimensionedVectorField>("g");
        return value.value();
    }
    else
    {
        return g_;
    }
}


void Foam::fv::castingSolidificationMelting::update
(
    const volScalarField& Cp
) const
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    if (debug)
    {
        Info<< indent
            << type() << ": " << name()
            << " - updating phase indicator" << endl;
    }

    // Update old-time alpha1 field once per physical time step.
    alpha1_.oldTime();

    const volScalarField& T = mesh().lookupObject<volScalarField>(TName_);

    const volScalarField* CarbonPtr = nullptr;

    // Gate 5A fallback field.
    // Used only when Carbon is not already owned by another model/function object.
    autoPtr<volScalarField> CarbonFilePtr;

    if (compositionDependentLiquidus_)
    {
        if (mesh().foundObject<volScalarField>(CarbonName_))
        {
            // Normal coupled case:
            // use the live transported Carbon field from the objectRegistry.
            CarbonPtr =
                &mesh().lookupObject<volScalarField>(CarbonName_);
        }
        else
        {
            // Standalone Gate 5A thermodynamic test:
            // read Carbon from the starting-time directory without registering it.
            const word CarbonInstance =
                mesh().time().timeName
                (
                    mesh().time().startTime().value()
                );

            CarbonFilePtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        CarbonName_,
                        CarbonInstance,
                        mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh()
                )
            );

            CarbonPtr = &CarbonFilePtr();
        }
    }

    zone_.regenerate();
    const labelList& cells = zone_.zone();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar Tc = T[celli];
        const scalar Cpc = Cp[celli];

        scalar TeqOld = 0;

        if (compositionDependentLiquidus_)
        {
            const scalar CarbonRaw = (*CarbonPtr)[celli];

            if (CarbonRaw < -SMALL)
            {
                FatalErrorInFunction
                    << "Negative Carbon concentration detected." << nl
                    << "    cell      = " << celli << nl
                    << "    Carbon    = " << CarbonRaw << nl
                    << "    fieldName = " << CarbonName_
                    << abort(FatalError);
            }

            const scalar Carbon = max(CarbonRaw, scalar(0));
            const scalar alphaOld = alpha1_[celli];

            const scalar denominatorOld =
                partitionCoefficient_
              + (1.0 - partitionCoefficient_)*alphaOld;

            const scalar CarbonLOld = Carbon/denominatorOld;

            TeqOld =
                Tliq_
              + liquidusSlope_*(CarbonLOld - CarbonRef_);
        }
        else
        {
            TeqOld =
                max
                (
                    Tsol_,
                    Tsol_
                  + (Tliq_ - Tsol_)
                   *(alpha1_[celli] - alpha1e_)
                   /(1 - alpha1e_)
                );
        }

        const scalar alpha1New =
            alpha1_[celli]
          + relax_*Cpc*(Tc - TeqOld)/L_;

        alpha1_[celli] = max(0, min(alpha1New, 1));

        if (compositionDependentLiquidus_)
        {
            const scalar Carbon = max((*CarbonPtr)[celli], scalar(0));
            const scalar alphaNew = alpha1_[celli];

            const scalar denominatorNew =
                partitionCoefficient_
              + (1.0 - partitionCoefficient_)*alphaNew;

            const scalar CarbonL = Carbon/denominatorNew;
            const scalar CarbonS = partitionCoefficient_*CarbonL;

            const scalar Teq =
                Tliq_
              + liquidusSlope_*(CarbonL - CarbonRef_);

            CarbonL_[celli] = CarbonL;
            CarbonS_[celli] = CarbonS;
            Teq_[celli] = Teq;
            deltaT_[i] = Tc - Teq;
        }
        else
        {
            const scalar Teq =
                max
                (
                    Tsol_,
                    Tsol_
                  + (Tliq_ - Tsol_)
                   *(alpha1_[celli] - alpha1e_)
                   /(1 - alpha1e_)
                );

            CarbonL_[celli] = 0;
            CarbonS_[celli] = 0;
            Teq_[celli] = Teq;
            deltaT_[i] = Tc - Teq;
        }
    }

    alpha1_.correctBoundaryConditions();
    CarbonL_.correctBoundaryConditions();
    CarbonS_.correctBoundaryConditions();
    Teq_.correctBoundaryConditions();

    curTimeIndex_ = mesh().time().timeIndex();
}


template<class RhoFieldType>
void Foam::fv::castingSolidificationMelting::apply
(
    const RhoFieldType& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< indent
            << type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    dimensionedScalar L("L", dimEnergy/dimMass, L_);

    const surfaceScalarField& phi =
        mesh().lookupObject<surfaceScalarField>(phiName_);

    // Conservative phase-change rate for a continuously moving strand.
    // The convective term was validated in Gate 2A/2B and is intentionally
    // retained unchanged for Gate 5A.
    if (eqn.psi().dimensions() == dimTemperature)
    {
        eqn -=
            L/Cp
           *(
                fvc::ddt(rho, alpha1_)
              + fvc::div(phi*fvc::interpolate(alpha1_))
            );
    }
    else
    {
        eqn -=
            L
           *(
                fvc::ddt(rho, alpha1_)
              + fvc::div(phi*fvc::interpolate(alpha1_))
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::castingSolidificationMelting::castingSolidificationMelting
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    Tsol_(NaN),
    Tliq_(NaN),
    alpha1e_(NaN),
    compositionDependentLiquidus_(false),
    CarbonName_("Carbon"),
    partitionCoefficient_(1.0),
    CarbonRef_(0.0),
    liquidusSlope_(0.0),
    L_(NaN),
    relax_(NaN),
    mode_(thermoMode::thermo),
    rhoRef_(NaN),
    TName_(word::null),
    CpName_(word::null),
    UName_(word::null),
    phiName_(word::null),
    Cu_(NaN),
    q_(NaN),
    pullVelocity_(vector::zero),
    beta_(NaN),
    alpha1_
    (
        IOobject
        (
            this->name() + ":alpha1",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
            this->name() + ":CarbonL",
            mesh.time().name(),
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
            this->name() + ":CarbonS",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Teq_
    (
        IOobject
        (
            this->name() + ":Teq",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    curTimeIndex_(-1),
    deltaT_(zone_.nCells(), 0)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::castingSolidificationMelting::addSupFields() const
{
    switch (mode_)
    {
        case thermoMode::thermo:
        {
            const basicThermo& thermo =
                mesh().lookupObject<basicThermo>(physicalProperties::typeName);

            return wordList({UName_, thermo.he().name()});
        }
        case thermoMode::lookup:
        {
            return wordList({UName_, TName_});
        }
    }

    return wordList::null();
}


void Foam::fv::castingSolidificationMelting::addSup
(
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    apply(geometricOneField(), eqn);
}


void Foam::fv::castingSolidificationMelting::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    apply(rho, eqn);
}


void Foam::fv::castingSolidificationMelting::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< indent
            << type() << ": applying source to " << eqn.psi().name() << endl;
    }

    const volScalarField Cp(this->Cp());

    update(Cp);

    const vector g = this->g();

    scalarField& Sp = eqn.diag();
    vectorField& Su = eqn.source();
    const scalarField& V = mesh().V();

    const labelList& cells = zone_.zone();

    forAll(cells, i)
    {
        const label celli = cells[i];

        const scalar Vc = V[celli];
        const scalar alpha1c = alpha1_[celli];

        const scalar S = -Cu_*sqr(1.0 - alpha1c)/(pow3(alpha1c) + q_);
        const vector Sb = rhoRef_*g*beta_*deltaT_[i];

        Sp[celli] += Vc*S;
        Su[celli] += Vc*(Sb + S*pullVelocity_);
    }
}


void Foam::fv::castingSolidificationMelting::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    addSup(U, eqn);
}


bool Foam::fv::castingSolidificationMelting::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::castingSolidificationMelting::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::castingSolidificationMelting::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::castingSolidificationMelting::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


bool Foam::fv::castingSolidificationMelting::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //