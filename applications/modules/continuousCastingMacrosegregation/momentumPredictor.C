/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    Derived from the OpenFOAM Foundation v14 incompressibleFluid solver
    module and distributed under the GNU General Public License.

\*---------------------------------------------------------------------------*/

#include "continuousCastingMacrosegregation.H"
#include "fvcGrad.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::continuousCastingMacrosegregation::momentumPredictor()
{
    volVectorField& U(U_);

    // Update thermo-solutal Boussinesq acceleration.
    updateBuoyancy();

    // Update moving-solid BKC inverse permeability and the corresponding
    // kinematic Darcy damping coefficient nu_l/K.
    updateBKCDrag(false);

    const dimensionedVector solidVelocityDim
    (
        "solidVelocity",
        dimLength/dimTime,
        solidVelocity_
    );

    // Moving-solid BKC drag:
    //
    //   a_D = -D (U-us),  D = nu_l/K
    //
    // The -D*U part is treated implicitly through Sp(D,U). The +D*us part
    // is an explicit source on the right-hand side. For solidVelocity=0 this
    // reduces exactly to the validated stationary-solid BKC formulation.
    tUEqn =
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + momentumTransport->divDevSigma(U)
      + fvm::Sp(bkcDragCoeff_, U)
     ==
        fvModels().source(U)
      + bkcDragCoeff_*solidVelocityDim
    );

    fvVectorMatrix& UEqn = tUEqn.ref();

    // For fvMatrix, operator-= with a VolField adds that explicit field to
    // the right-hand-side source. buoyancyTotal_ is acceleration because
    // this solver uses kinematic pressure.
    UEqn -= buoyancyTotal_;

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints().constrain(U);
    }
}


// ************************************************************************* //
