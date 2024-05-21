/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

/*---------------------------------------------------------------------------* \
  Stright copy of kOmegaSST with all important operations in one file.
\*---------------------------------------------------------------------------*/

#include "frozenkOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenkOmegaSST<BasicTurbulenceModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_LES_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_LES_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenkOmegaSST<BasicTurbulenceModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_LES_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenkOmegaSST<BasicTurbulenceModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
frozenkOmegaSST<BasicTurbulenceModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void frozenkOmegaSST<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_LES_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void frozenkOmegaSST<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
frozenkOmegaSST<BasicTurbulenceModel>::Pk
(
    const volScalarField::Internal& G
) const
{
  return min(G, (c1_*betaStar_)*k_LES_()*omega_());
}


template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
frozenkOmegaSST<BasicTurbulenceModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
frozenkOmegaSST<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_LES_,
            dimVolume*this->rho_.dimensions()*k_LES_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
frozenkOmegaSST<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> frozenkOmegaSST<BasicTurbulenceModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
frozenkOmegaSST<BasicTurbulenceModel>::frozenkOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    //========================== Frozen parameters ======================
    // Set these to 0 to turn-off influence of LES, 1 to include fully.
    // Note: LES velocity field and k are always used in the omega equation.
    useRST_   // Use the LES RST (via a_ij) in the omega-eqn production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "useRST",
            this->coeffDict_,
            1.0
        )
    ),
    usekDeficit_  // Use the k-eqn residual to correct the omega production
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "usekDeficit",
            this->coeffDict_,
            1.0
        )
    ),
    // Ramping gradually introduce corrections (both kDeficit and bDelta) to
    // aid solver stability.  Before `rampStartTime` corrections are zero,
    // after `rampEndTime` they are 1.0.  Linear in between.  Default is
    // full correction from beginning.
    rampStartTime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "rampStartTime",
            this->coeffDict_,
	    dimTime,
            0
        )
    ),
    rampEndTime_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "rampEndTime",
            this->coeffDict_,
	    dimTime,
            1
        )
    ),
    xi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "xi_ramp",
            this->coeffDict_,
	    dimless,
            1
        )
    ),
    //========================== LES fields ==============================
    k_LES_
    (
        IOobject
        (
            IOobject::groupName("k_LES", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    tauij_LES_
    (
        IOobject
        (
            "tauij_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    aij_LES_
    (
        IOobject
        (
            "aij_LES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tauij_LES_ - ((2.0/3.0)*I)*k_LES_
    ),

    //========================== Unknown fields ============================
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kDeficit_
    (
        IOobject(
	    "kDeficit",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("kDeficit", dimensionSet(0,2,-3,0,0,0,0), 0.0)
    ),
    aijDelta_
    (
        IOobject
        (
            "aijDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*symm(fvc::grad(this->U_))*this->nut_
    ),
    bijDelta_
    (
        IOobject
        (
            "bijDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*symm(fvc::grad(this->U_))/omega_
    ),

     
    
    //========================== Misc fields ============================
    aijBoussinesq_
    (
        IOobject
        (
            "aijBoussinesq",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*symm(fvc::grad(this->U_))*this->nut_
    ),
    Pk_
    (
        IOobject
        (
            "Pk",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Pk", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),
    PkBoussinesq_
    (
        IOobject
        (
            "PkBoussinesq",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("PkBoussinesq", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),
    PkDelta_
    (
        IOobject
        (
            "PkDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("PkDelta", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),
    PkLES_
    (
        IOobject
        (
            "PkLES",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("PkLES", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    ),
    y_(
       IOobject
       (
            "walldist",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       wallDist::New(this->mesh_).y()
    ),
    gradU_
    (
        IOobject
        (
            "gradU",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(this->U_)
    ),
    gradk_
    (
        IOobject
        (
            "gradk",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(this->k_LES_)  // Only compute once - k=k_LES is fixed in frozen
    ),
     
    gradomega_
    (
       IOobject 
       (
           "gradomega",
           this->runTime_.timeName(),
           this->mesh_,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       this->mesh_,
       dimensionedVector("gradomega", dimensionSet(0,-1,-1,0,0,0,0), Zero)
    ),
    nu_
    (
        IOobject
        (
            "nuOut",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->nu()//.internalField()
    )   

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
    bound(k_LES_, this->kMin_);
    bound(omega_, this->omegaMin_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool frozenkOmegaSST<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix> 
frozenkOmegaSST<BasicTurbulenceModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: frozenkOmegaSST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)    // linear part
      + fvc::div(dev(this->aij_LES_))
    );
}
 

//- Return the modified effective stress tensor  
template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
frozenkOmegaSST<BasicTurbulenceModel>::devRhoReff() const
{
    Info << "In: frozenkOmegaSST::devRhoReff()" << endl;
    return volSymmTensorField::New
    (
        IOobject::groupName("devRhoReff", this->alphaRhoPhi_.group()),
        (-(this->alpha_*this->rho_*this->nu()))
       *dev(twoSymm(fvc::grad(this->U_)))
       + dev(this->aij_LES_) 
    );
}

//- Return the modified source term for the momentum equation
template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
frozenkOmegaSST<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    Info << "In: frozenkOmegaSST::divDevRhoReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nu())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nu(), U)
      + fvc::div(dev(this->aij_LES_))
    );
}

//- Return the modified source term for the momentum equation
template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
frozenkOmegaSST<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    Info << "In: frozenkOmegaSST::divDevRhoReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*rho*this->nu())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*rho*this->nu(), U)
      + fvc::div(dev(aij_LES_))
    );
}  


  
template<class BasicTurbulenceModel>
void frozenkOmegaSST<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_) {    return;   }
    Info << "In: frozenkOmegaSST.correct()" << endl;
    
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    volScalarField& omega_ = this->omega_;
    volScalarField& k_LES_ = this->k_LES_;

    BasicTurbulenceModel::correct();

    // This is the current iteration, 1000, etc.
    const dimensionedScalar time = this->runTime_;
    xi_ = (time < rampStartTime_)? 0.0:
          (time > rampEndTime_)? 1.0:
          (time - rampStartTime_) / (rampEndTime_ - rampStartTime_);
    Info << "Corrections: xi = " << xi_.value() <<
          ", kDeficit factor = " << (xi_*usekDeficit_).value() <<
          ", bijDelta factor = " << (xi_*useRST_).value() << endl;

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );


    tmp<volTensorField> tgradU = fvc::grad(U);
    volSymmTensorField S(symm(tgradU()) / omega_);
    volTensorField W(skew(tgradU()) / omega_);
    volScalarField S2(2*magSqr(symm(tgradU())));

    volScalarField GbyNu(dev(twoSymm(tgradU())) && tgradU());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    volScalarField G2
    (
         "G2",
	 nut*GbyNu - useRST_ * xi_ * (aijDelta_ && tgradU())
    );

    // Write all these production terms to output
    PkLES_        = -tauij_LES_ && tgradU();
    Pk_           = G2;
    PkBoussinesq_ = nut*GbyNu;
    PkDelta_      = -aijDelta_ && tgradU();

    // "Free" temporary variable
    tgradU.clear();

    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k_LES_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent specific dissipation rate equation - compared with SSTm from:
        //   https://turbmodels.larc.nasa.gov/sst.html
        tmp<fvScalarMatrix> omegaEqn
        (
	 // 1. Time-derivative = d/dt(rho omega)
	 fvm::ddt(alpha, rho, omega_) 
	 // 2. Advection of omaga = d/dx_j(rho u_j omega)
	 + fvm::div(alphaRhoPhi, omega_)
	 // 3. Molecular + turb. diffusion
	 //    Where DomegaEff(F1) = alphaOmega(F1)*this->nut_ + this->nu() (kOmegaSSTBase.H)
	 - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega_)
         ==
	 // 4. Production of omega
	 //    This term employs the production limiter min(P, 20 beta* rho omega k), which
	 //    when using the definition of mu_t, gives the following min-max relation.
	 //    I.e. consistent with NASA SSTm.
	 //
	 //    In Martin's original code the correction to the production due to bijDelta
	 //    was included in the limiting, but the kDefecit correction not included
	 //    (unlimited).  We follow this pattern here.
	   alpha()*rho()*gamma
           *min
            (
	     G2 / nut(),
                (this->c1_/this->a1_)*this->betaStar_*omega_()
               *max(this->a1_*omega_(), this->b1_*F23()*sqrt(S2()))
            )
	 // 4b. kDeficit term (error in k eqn transferred to omega)
	 //+ alpha()*rho()*gamma * kDeficit_/nut()*(usekDeficit_*xi_) //Removed due to instability
	 // 5. Compressible part of production (for divU != 0)
	 - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
	 // 6. Dissipation of omega = -beta rho omega^2
	 - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
	 // 7. Cross-diffusion = +2(1-F1)...
	 //    Where CDkOmega() = (2*alphaOmega2)*(grad(k) grad(omega))/omega
	 - fvm::SuSp
	 (
	  alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
	  omega_
	  )
	 // 8. Not clear - Scale-adaptive simulation term?  Inactive = 0 here?
	 //    Defined in (kOmegaSSTBase.C)
	 + this->Qsas(S2(), gamma, beta)
	 // 9. Case-specific omega sources = 0 (check?)
	 + this->omegaSource()
	 // 10. Modification defined in parameter files = 0
	 + fvOptions(alpha, rho, omega_)
        );
	omega_.boundaryFieldRef().updateCoeffs();

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);

	// Update gradomega (for output only)
	gradomega_ = fvc::grad(omega_);
    }

    /// // Turbulent kinetic energy equation
    /// tmp<fvScalarMatrix> kEqn
    /// (
    ///     fvm::ddt(alpha, rho, k_)
    ///   + fvm::div(alphaRhoPhi, k_)
    ///   - fvm::laplacian(alpha*rho*this->DkEff(F1), k_)
    ///  ==
    ///     alpha()*rho()*G2  //this->Pk(G)
    ///   - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
    ///   - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, F23), k_)
    ///   + this->kSource()
    ///   + fvOptions(alpha, rho, k_)
    /// );
    /// 
    /// kEqn.ref().relax();
    /// fvOptions.constrain(kEqn.ref());
    /// solve(kEqn);
    /// fvOptions.correct(k_);
    /// bound(k_LES_, this->kMin_);
    
    // Notes:
    //   - Incompressible only due to the change from fvm to fvc (why?)
    //   - Is ddt necessary here in the steady case?
    //   - kSource() and fvOptions() not included
    //   - .ref() returns a non-const reference to the internalField, whereas
    //     both v() and .internalField(), as well as () return const references.
    kDeficit_.ref() = //  fvc::ddt(alpha, rho, k_LES_)
      fvc::div(alphaRhoPhi, k_LES_)()()
      - fvc::laplacian(alpha*rho*this->DkEff(F1), k_LES_)()()
      - alpha()*rho()*this->Pk(G2())
      + (2.0/3.0)*alpha()*rho()*divU*k_LES_()
      + alpha()*rho()*this->epsilonByk(F1, F23) * k_LES_();
   
    // Re-calculate viscosity
    this->correctNut(S2, F23);

    // Re-calculate aijDelta
    aijBoussinesq_  = -this->nut_ * twoSymm(fvc::grad(U));  // Boussinesq approx of <u_i' u_j'>
    aijDelta_       = aij_LES_ - aijBoussinesq_;
    bijDelta_.ref() = aijDelta_() / (2*(k_LES_() + this->kMin_));
    //Info << "RPD" << endl;
}



// /////////////////// Original SST correct() ///////////////////////
// 
// template<class BasicTurbulenceModel>
// void frozenkOmegaSST<BasicTurbulenceModel>::correct()
// {
//     if (!this->turbulence_)
//     {
//         return;
//     }
// 
//     // Local references
//     const alphaField& alpha = this->alpha_;
//     const rhoField& rho = this->rho_;
//     const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
//     const volVectorField& U = this->U_;
//     volScalarField& nut = this->nut_;
//     fv::options& fvOptions(fv::options::New(this->mesh_));
// 
//     BasicTurbulenceModel::correct();
// 
//     volScalarField::Internal divU
//     (
//         fvc::div(fvc::absolute(this->phi(), U))()()
//     );
// 
//     tmp<volTensorField> tgradU = fvc::grad(U);
//     volScalarField S2(2*magSqr(symm(tgradU())));
//     volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
//     volScalarField::Internal G(this->GName(), nut()*GbyNu);
//     tgradU.clear();
// 
//     // Update omega and G at the wall
//     omega_.boundaryFieldRef().updateCoeffs();
// 
//     volScalarField CDkOmega
//     (
//         (2*alphaOmega2_)*(fvc::grad(k_LES_) & fvc::grad(omega_))/omega_
//     );
// 
//     volScalarField F1(this->F1(CDkOmega));
//     volScalarField F23(this->F23());
// 
//     {
//         volScalarField::Internal gamma(this->gamma(F1));
//         volScalarField::Internal beta(this->beta(F1));
// 
//         // Turbulent frequency equation
//         tmp<fvScalarMatrix> omegaEqn
//         (
//             fvm::ddt(alpha, rho, omega_)
//           + fvm::div(alphaRhoPhi, omega_)
//           - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
//          ==
//             alpha()*rho()*gamma
//            *min
//             (
//                 GbyNu,
//                 (c1_/a1_)*betaStar_*omega_()
//                *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
//             )
//           - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
//           - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
//           - fvm::SuSp
//             (
//                 alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
//                 omega_
//             )
//           + Qsas(S2(), gamma, beta)
//           + omegaSource()
//           + fvOptions(alpha, rho, omega_)
//         );
// 
//         omegaEqn.ref().relax();
//         fvOptions.constrain(omegaEqn.ref());
//         omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
//         solve(omegaEqn);
//         fvOptions.correct(omega_);
//         bound(omega_, this->omegaMin_);
//     }
// 
//     // Turbulent kinetic energy equation
//     tmp<fvScalarMatrix> kEqn
//     (
//         fvm::ddt(alpha, rho, k_LES_)
//       + fvm::div(alphaRhoPhi, k_LES_)
//       - fvm::laplacian(alpha*rho*DkEff(F1), k_LES_)
//      ==
//         alpha()*rho()*Pk(G)
//       - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_LES_)
//       - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_LES_)
//       + kSource()
//       + fvOptions(alpha, rho, k_LES_)
//     );
// 
//     kEqn.ref().relax();
//     fvOptions.constrain(kEqn.ref());
//     solve(kEqn);
//     fvOptions.correct(k_LES_);
//     bound(k_LES_, this->kMin_);
// 
//     correctNut(S2, F23);
// }
// 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
