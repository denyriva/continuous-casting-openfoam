/*---------------------------------------------------------------------------*\
  DongLaunderSharmaKE v3 diagnostic
  OpenFOAM Foundation v14-compatible custom RAS model.

  Base:
      Foam::RASModels::LaunderSharmaKE

  Dong-style solidification damping:
      a_turb = nu * (180/lambda2^2)
               * (1 - fl^2)/(fl^3 + xi)

      S_k       = sign * scale * a_turb * k
      S_epsilon = sign * scale * a_turb * epsilon

  Default sign = -1 => damping/sink.
\*---------------------------------------------------------------------------*/

#ifndef DongLaunderSharmaKE_H
#define DongLaunderSharmaKE_H

#include "LaunderSharmaKE.H"
#include "fvmSup.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class DongLaunderSharmaKE
:
    public LaunderSharmaKE<BasicMomentumTransportModel>
{
    typedef LaunderSharmaKE<BasicMomentumTransportModel> Base;

    scalar solidificationSourceSign_;
    scalar solidificationSourceScale_;
    scalar lambda2_;
    scalar solidificationXi_;
    word solidFractionFieldName_;

    // Cell-centred consistency diagnostic for the inherited
    // Launder-Sharma nut update.  This does not alter the model.
    bool writeNutConsistencyDiagnostics_;
    scalar nutDiagnosticFloor_;

    tmp<volScalarField> solidificationDampingRate() const
    {
        if (!this->mesh_.template foundObject<volScalarField>
            (solidFractionFieldName_))
        {
            FatalErrorInFunction
                << "Required solid-fraction field '"
                << solidFractionFieldName_
                << "' is not present in the object registry."
                << exit(FatalError);
        }

        const volScalarField& fs =
            this->mesh_.template lookupObject<volScalarField>
            (solidFractionFieldName_);

        const volScalarField fl
        (
            volScalarField::New
            (
                "DongLiquidFraction",
                min(max(scalar(1) - fs, scalar(0)), scalar(1))
            )
        );

        const dimensionedScalar Am
        (
            "Am",
            dimless/sqr(dimLength),
            180.0/sqr(lambda2_)
        );

        return
            this->nu()
           *Am
           *(scalar(1) - sqr(fl))
           /(pow3(fl) + solidificationXi_);
    }


protected:

    virtual tmp<fvScalarMatrix> kSource() const
    {
        tmp<fvScalarMatrix> tSource = Base::kSource();

        if (mag(solidificationSourceSign_) <= SMALL
         || solidificationSourceScale_ <= SMALL)
        {
            return tSource;
        }

        const tmp<volScalarField> tRate = solidificationDampingRate();

        if (solidificationSourceSign_ < 0)
        {
            tSource.ref() -=
                fvm::Sp
                (
                    (-solidificationSourceSign_)
                   *solidificationSourceScale_
                   *tRate(),
                    this->k_
                );
        }
        else
        {
            tSource.ref() +=
                fvm::Sp
                (
                    solidificationSourceSign_
                   *solidificationSourceScale_
                   *tRate(),
                    this->k_
                );
        }

        return tSource;
    }


    virtual tmp<fvScalarMatrix> epsilonSource() const
    {
        tmp<fvScalarMatrix> tSource = Base::epsilonSource();

        if (mag(solidificationSourceSign_) <= SMALL
         || solidificationSourceScale_ <= SMALL)
        {
            return tSource;
        }

        const tmp<volScalarField> tRate = solidificationDampingRate();

        if (solidificationSourceSign_ < 0)
        {
            tSource.ref() -=
                fvm::Sp
                (
                    (-solidificationSourceSign_)
                   *solidificationSourceScale_
                   *tRate(),
                    this->epsilon_
                );
        }
        else
        {
            tSource.ref() +=
                fvm::Sp
                (
                    solidificationSourceSign_
                   *solidificationSourceScale_
                   *tRate(),
                    this->epsilon_
                );
        }

        return tSource;
    }


public:

    TypeName("DongLaunderSharmaKE");

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    DongLaunderSharmaKE
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    )
    :
        Base(alpha, rho, U, alphaRhoPhi, phi, viscosity, type),

        solidificationSourceSign_
        (
            this->typeDict(type).template lookupOrDefault<scalar>
            (
                "solidificationSourceSign",
               -1.0
            )
        ),

        solidificationSourceScale_
        (
            this->typeDict(type).template lookupOrDefault<scalar>
            (
                "solidificationSourceScale",
                1.0
            )
        ),

        lambda2_
        (
            this->typeDict(type).template lookupOrDefault<scalar>
            (
                "lambda2",
                1e-4
            )
        ),

        solidificationXi_
        (
            this->typeDict(type).template lookupOrDefault<scalar>
            (
                "solidificationXi",
                1e-3
            )
        ),

        solidFractionFieldName_
        (
            this->typeDict(type).template lookupOrDefault<word>
            (
                "solidFractionField",
                "fs"
            )
        ),

        writeNutConsistencyDiagnostics_
        (
            this->typeDict(type).template lookupOrDefault<bool>
            (
                "writeNutConsistencyDiagnostics",
                false
            )
        ),

        nutDiagnosticFloor_
        (
            this->typeDict(type).template lookupOrDefault<scalar>
            (
                "nutDiagnosticFloor",
                1e-30
            )
        )
    {
        if (lambda2_ <= SMALL)
        {
            FatalErrorInFunction
                << "lambda2 must be > 0. Current value = "
                << lambda2_ << " m"
                << exit(FatalError);
        }

        if (solidificationXi_ <= 0)
        {
            FatalErrorInFunction
                << "solidificationXi must be > 0. Current value = "
                << solidificationXi_
                << exit(FatalError);
        }

        Info<< "DongLaunderSharmaKE v2:" << nl
            << "    base model               = LaunderSharmaKE" << nl
            << "    solidification sign      = "
            << solidificationSourceSign_ << nl
            << "    solidification scale     = "
            << solidificationSourceScale_ << nl
            << "    lambda2                  = "
            << lambda2_ << " m" << nl
            << "    solidificationXi         = "
            << solidificationXi_ << nl
            << "    solidFractionField       = "
            << solidFractionFieldName_ << nl
            << "    write nut diagnostics    = "
            << writeNutConsistencyDiagnostics_ << nl
            << "    nut diagnostic floor     = "
            << nutDiagnosticFloor_ << " m2/s" << nl
            << "    damping form             = "
            << "nu*(180/lambda2^2)*(1-fl^2)/(fl^3+xi)"
            << endl;
    }



    virtual void correct()
    {
        // Run the unmodified Foundation Launder-Sharma update first:
        // epsilon solve -> epsilon bound -> k solve -> k bound -> correctNut().
        Base::correct();

        if
        (
            !writeNutConsistencyDiagnostics_
         || !this->runTime_.writeTime()
        )
        {
            return;
        }

        const volScalarField fMuExpected
        (
            IOobject
            (
                "fMuExpected",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            exp
            (
               -3.4
               /sqr
                (
                    scalar(1)
                  + sqr(this->k_)
                   /(this->nu()*this->epsilon_)/50.0
                )
            )
        );

        const volScalarField nutExpected
        (
            IOobject
            (
                "nutExpected",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            0.09*fMuExpected*sqr(this->k_)/this->epsilon_
        );

        const dimensionedScalar nutFloor
        (
            "nutDiagnosticFloor",
            this->nut_.dimensions(),
            nutDiagnosticFloor_
        );

        const volScalarField nutDifference
        (
            IOobject
            (
                "nutDifference",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->nut_ - nutExpected
        );

        const volScalarField nutRatio
        (
            IOobject
            (
                "nutRatio",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->nut_/max(nutExpected, nutFloor)
        );

        const volScalarField nutRelativeMismatch
        (
            IOobject
            (
                "nutRelativeMismatch",
                this->runTime_.name(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mag(nutDifference)
           /max
            (
                max(mag(this->nut_), mag(nutExpected)),
                nutFloor
            )
        );

        fMuExpected.write();
        nutExpected.write();
        nutDifference.write();
        nutRatio.write();
        nutRelativeMismatch.write();

        Info<< "DongLaunderSharmaKE nut consistency diagnostic" << nl
            << "    time                  = "
            << this->runTime_.value() << " s" << nl
            << "    max(nut actual)       = "
            << gMax(this->nut_.primitiveField()) << " m2/s" << nl
            << "    max(nut expected)     = "
            << gMax(nutExpected.primitiveField()) << " m2/s" << nl
            << "    max |nut-expected|    = "
            << gMax(mag(nutDifference.primitiveField())) << " m2/s" << nl
            << "    max relative mismatch = "
            << gMax(nutRelativeMismatch.primitiveField()) << nl
            << "    wrote fields: fMuExpected nutExpected nutDifference "
            << "nutRatio nutRelativeMismatch"
            << endl;
    }


    virtual bool read()
    {
        if (!Base::read())
        {
            return false;
        }

        const dictionary& coeffs = this->typeDict();

        solidificationSourceSign_ =
            coeffs.template lookupOrDefault<scalar>
            (
                "solidificationSourceSign",
                solidificationSourceSign_
            );

        solidificationSourceScale_ =
            coeffs.template lookupOrDefault<scalar>
            (
                "solidificationSourceScale",
                solidificationSourceScale_
            );

        lambda2_ =
            coeffs.template lookupOrDefault<scalar>
            (
                "lambda2",
                lambda2_
            );

        solidificationXi_ =
            coeffs.template lookupOrDefault<scalar>
            (
                "solidificationXi",
                solidificationXi_
            );

        solidFractionFieldName_ =
            coeffs.template lookupOrDefault<word>
            (
                "solidFractionField",
                solidFractionFieldName_
            );

        writeNutConsistencyDiagnostics_ =
            coeffs.template lookupOrDefault<bool>
            (
                "writeNutConsistencyDiagnostics",
                writeNutConsistencyDiagnostics_
            );

        nutDiagnosticFloor_ =
            coeffs.template lookupOrDefault<scalar>
            (
                "nutDiagnosticFloor",
                nutDiagnosticFloor_
            );

        return true;
    }
};

} // End namespace RASModels
} // End namespace Foam

#endif
