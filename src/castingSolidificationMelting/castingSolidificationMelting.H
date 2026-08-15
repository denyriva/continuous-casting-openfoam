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

Class
    Foam::fv::castingSolidificationMelting

Description
    Continuous-casting solidification/melting fvModel derived from the
    Foundation solidificationMelting model.

    The model includes:

      - latent heat of fusion;
      - mushy-zone momentum resistance;
      - prescribed pull velocity for the solidifying strand;
      - thermal buoyancy;
      - conservative transport of latent energy;
      - optional composition-dependent binary phase equilibrium.

    With compositionDependentLiquidus = false, the original fixed-Tsol/Tliq
    phase-fraction closure is retained.

    With compositionDependentLiquidus = true, the model reads the conserved
    mixture carbon field Carbon and evaluates, locally,

        CarbonL = Carbon/[kp + (1-kp)*alpha1]
        CarbonS = kp*CarbonL
        Teq     = TliqRef + liquidusSlope*(CarbonL - CarbonRef)

    where TliqRef is the dictionary entry Tliq.  The phase fraction is updated
    with the same enthalpy-correction algorithm used by the original model,
    but using this composition-dependent equilibrium temperature.

    The composition-dependent closure implemented here is a linearized binary
    equilibrium model with constant partition coefficient.  It does not add
    transport equations for alpha1, CarbonL, CarbonS, or Teq.  Carbon remains
    the transported/conserved solute variable.

    Diagnostic fields written by the model are:

        <modelName>:alpha1
        <modelName>:CarbonL
        <modelName>:CarbonS
        <modelName>:Teq

Usage
    Legacy mode:

    \verbatim
    solidificationMelting1
    {
        type            castingSolidificationMelting;
        cellZone        all;

        Tsol            1593.15;
        Tliq            1688.15;
        alpha1e         0;

        L               220000;
        relax           0.9;

        thermoMode      thermo;
        rhoRef          7035;
        beta            2.4e-4;

        Cu              1e5;
        q               1e-3;
        pullVelocity    (0 0 0.013333333);
    }
    \endverbatim

    Gate 5A composition-dependent mode:

    \verbatim
    solidificationMelting1
    {
        type            castingSolidificationMelting;
        cellZone        all;

        Tsol            1593.15;   // Used by legacy mode
        Tliq            1688.15;   // Reference liquidus at CarbonRef
        alpha1e         0;

        compositionDependentLiquidus true;
        Carbon          Carbon;
        partitionCoefficient 0.25;
        CarbonRef       0.01;
        liquidusSlope   -1000;     // K per unit mass fraction; test value

        L               220000;
        relax           0.9;

        thermoMode      thermo;
        rhoRef          7035;
        beta            0;

        Cu              1e5;
        q               1e-3;
        pullVelocity    (0 0 0);
        g               (0 0 0);
    }
    \endverbatim

SourceFiles
    castingSolidificationMelting.C

\*---------------------------------------------------------------------------*/

#ifndef castingSolidificationMelting_H
#define castingSolidificationMelting_H

#include "fvModel.H"
#include "fvCellZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class castingSolidificationMelting Declaration
\*---------------------------------------------------------------------------*/

class castingSolidificationMelting
:
    public fvModel
{
public:

    enum class thermoMode
    {
        thermo,
        lookup
    };

    static const NamedEnum<thermoMode, 2> thermoModeTypeNames_;


private:

    // Private Data

        //- The cellZone the fvModel applies to
        mutable fvCellZone zone_;

        //- Solidus temperature used by the legacy fixed phase diagram [K]
        scalar Tsol_;

        //- Legacy liquidus temperature, or reference liquidus temperature
        //  at CarbonRef_ when compositionDependentLiquidus_ is enabled [K]
        scalar Tliq_;

        //- Max eutectic melt fraction for the legacy closure [0-1[
        scalar alpha1e_;

        //- Enable the Gate 5A composition-dependent binary closure
        bool compositionDependentLiquidus_;

        //- Name of conserved mixture carbon mass-fraction field
        word CarbonName_;

        //- Equilibrium partition coefficient kp = CarbonS/CarbonL
        scalar partitionCoefficient_;

        //- Reference liquid carbon mass fraction associated with Tliq_
        scalar CarbonRef_;

        //- Liquidus slope [K per unit mass fraction]
        scalar liquidusSlope_;

        //- Latent heat of fusion [J/kg]
        scalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Thermodynamics mode
        thermoMode mode_;

        //- Reference density - typically the solid density
        scalar rhoRef_;

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Name of specific heat capacity field - default = "Cp" (optional)
        word CpName_;

        //- Name of velocity field - default = "U" (optional)
        word UName_;

        //- Name of flux field - default = "phi" (optional)
        word phiName_;

        //- Mushy region momentum sink coefficient [1/s]; default = 10^5
        scalar Cu_;

        //- Coefficient used in porosity calc - default = 0.001
        scalar q_;

        //- Velocity imposed on the solid phase [m/s]
        vector pullVelocity_;

        //- Reference temperature for thermal buoyancy [K]
        scalar TRef_;

        //- Thermal expansion coefficient [1/K]
        scalar beta_;

        //- Solutal expansion coefficient [per unit carbon mass fraction]
        scalar betaC_;

        //- Reference heat capacity for thermoMode::lookup
        scalar CpRef_;

        //- Gravitational acceleration
        vector g_;

        //- Phase fraction indicator field
        mutable volScalarField alpha1_;

        //- Gate 5A diagnostic: equilibrium liquid carbon concentration
        mutable volScalarField CarbonL_;

        //- Gate 5A diagnostic: equilibrium solid carbon concentration
        mutable volScalarField CarbonS_;

        //- Gate 5A diagnostic: local equilibrium temperature [K]
        mutable volScalarField Teq_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;

        //- Temperature change cached for source calculation when alpha1 updated
        mutable scalarField deltaT_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the specific heat capacity field
        tmp<volScalarField> Cp() const;

        //- Return the gravity vector
        vector g() const;

        //- Update the model
        void update(const volScalarField& Cp) const;

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        inline void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn) const;


public:

    //- Runtime type information
    TypeName("castingSolidificationMelting");


    // Constructors

        //- Construct from explicit source name and mesh
        castingSolidificationMelting
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        castingSolidificationMelting
        (
            const castingSolidificationMelting&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluation

            //- Add explicit contribution to enthalpy equation
            virtual void addSup
            (
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add explicit contribution to compressible enthalpy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const castingSolidificationMelting&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //