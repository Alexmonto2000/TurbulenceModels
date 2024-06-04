/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "kOmegaSSTML_SpaRTA_V2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::kOmegaSSTML_SpaRTA_V2
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
    kOmegaSSTBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
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
    
    
    //====================== Propagation parameters ======================
    // Set these to 0 to turn off influence of frozen corrections, 1 to include fully.
    // The kDeficit and bDelta fields must come from the frozen solver.  Values in the
    // interval [0,1] make sense too.
    
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

    //========================== Fields from frozen solver ================
    epsilon_
    (
        IOobject(
	    "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("epsilon", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
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
        this->mesh_
    ),
    bijDelta_
    (
        IOobject
        (
            "bijDelta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
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
       this->mesh_,
       dimensionedTensor("gradU", dimensionSet(0,0,-1,0,0,0,0), Zero)
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
       this->mesh_,
       dimensionedVector("gradk", dimensionSet(0,1,-2,0,0,0,0), Zero)
    ),
    gradp_
    (
        IOobject
        (
            "gradp",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
       this->mesh_,
       dimensionedVector("gradp", dimensionSet(0,1,-2,0,0,0,0), Zero)
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

    Dk_
    (
        IOobject
        (
            "Dk",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Dk", dimensionSet(0,2,-3,0,0,0,0), 0.0) 
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

template<class BasicTurbulenceModel>
tmp<fvVectorMatrix> kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA_V2::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U)))) // linear part
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part 
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)  // non-linear part
    );
    //return
    //(
    //  - div(nuEff * grad(U)^T) // linear part (explicit fvc)
    //  - div(nuEff * grad(U))   // linear part (implicit fvm)
    //  + div(2.*k*bijDelta)     // non-linear part
    //);
}

//- Return the modified effective stress tensor  
template<class BasicTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::devRhoReff() const
{
    Info << "In: kOmegaSSTML_SpaRTA_V2::devRhoReff()" << endl;
    return volSymmTensorField::New
    (
        IOobject::groupName("devRhoReff", this->alphaRhoPhi_.group()),
        (-(this->alpha_*this->rho_*this->nuEff())) // 
       *dev(twoSymm(fvc::grad(this->U_)))
       + dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_
    );
}

//- Return the modified source term for the momentum equation
template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA_V2::divDevRhoReff(U)" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U)))) // 
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U) // 
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)
    );
}

//- Return the modified source term for the momentum equation
template<class BasicTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA_V2::divDevRhoReff(rho, U)" << endl;
    return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U)))) //
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U) //
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)
    );
}  
  
template<class BasicTurbulenceModel>
void kOmegaSSTML_SpaRTA_V2<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    Info << "In: kOmegaSSTML_SpaRTA_V2.correct()" << endl;
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    //Additional dictionary for magnetic field definition
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            this->runTime_.constant(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    //Kinematic viscosity
    dimensionedScalar nu
    (
        "nu",
        dimensionSet(0, 2,-1, 0, 0,0,0),
        transportProperties
    );
    //Density
    dimensionedScalar rhol
    (
        "rhol",
        dimensionSet(1, -3, 0, 0, 0,0,0),
        transportProperties
    );

    // Update omega and G at the wall
    volScalarField& omega_ = this->omega_;
    volScalarField& k_ = this->k_;

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

    //Custom code to ensure the most updated gradient of p is outputted
    //TODO: Now calculates gradp each iteration; can be optimized to only
    //calculate it each write time.
    word pName_ = ("p");
    tmp<volScalarField> p_ = U.db().lookupObject<volScalarField>(pName_);    
    gradp_ = fvc::grad(p_);  

    // Same but for Lorentz force
    // word lName_ = ("lorentz");
    // tmp<volVectorField> lorentz = U.db().lookupObject<volVectorField>(lName_);    
        
    tmp<volTensorField> tgradU = fvc::grad(U);
    volSymmTensorField S(symm(tgradU()) / omega_);
    volTensorField W(skew(tgradU()) / omega_);
    volScalarField S2(2*magSqr(symm(tgradU()))); 
    /*******************************Feature Addition**************************/
    //add features for models
    // volSymmTensorField T1(S);
    // volSymmTensorField T2((symm((S & W) - (W & S))));
    // volSymmTensorField T3((symm((S & S) - 0.3333*tr(S & S)*I)));
    // volSymmTensorField T4((symm((W & W) - 0.3333*tr(W & W)*I)));
    gradk_ = fvc::grad(k_);
    
    
    
    epsilon_ = (omega_*k_);
    // volScalarField q_gamma((mag(tgradU())*k_) / max(epsilon_, dimensionedScalar("small", dimensionSet(0,2,-3,0,0,0,0), 1.0e-10)  )); //ask about this this->k_
    // volScalarField q_Re( 2 - min(sqrt(k_)*y_ / (50 * nu),dimensionedScalar("small", dimensionSet(0,0,0,0,0,0,0), 2)));
    volScalarField Ak2( -2 * (k_ / pow(max(epsilon_, dimensionedScalar("small", dimensionSet(0,2,-3,0,0,0,0), 1.0e-5)) ,2)) * magSqr(gradk_));
    volScalarField q_T(mag(S));
    // volScalarField q_AlfS( (pow(k_ / epsilon_ , 1.5) / (rhol* sqrt(nu)))* 1.41421356 * mag(lorentz())/max(mag(S), dimensionedScalar("small", dimensionSet(0,0,0,0,0,0,0), 1.0e-10)  ));
/************************************************************************/


/*******************************Model Addition**************************/
    //Model kDeficit_
    // kDeficit_ = ((0.33667425 * (q_Re/0.76262944)*(epsilon_/0.072954023) +
    //  0.15326155 * sqrt(q_T/0.26157319)*(epsilon_/0.072954023) +
    //  0.84000332 * sqrt(q_AlfS/7.6895703)*(epsilon_/0.072954023)) * 0.0048252824); //ask about constant = 1
    kDeficit_ = ((0.54890036 *(epsilon_ /0.13471841) + 
    0.32530213 * exp((Ak2/0.017511644))*(epsilon_ /0.13471841) + 
    0.44071346 * ((q_T/0.09867829)/ (1 + pow(q_T/0.09867829, 2)))*(epsilon_ /0.13471841)) * 0.010597859);

    kDeficit_ = max(kDeficit_ , dimensionedScalar("small", dimensionSet(0,2,-3,0,0,0,0), 1.0e-8));
/************************************************************************/

    // lorentz.clear();

    volScalarField GbyNu(dev(twoSymm(tgradU())) && tgradU());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    volScalarField G2
    (
         "G2",
     nut*GbyNu - xi_ * useRST_ * (2*(this->k_)*bijDelta_ && tgradU()) // 
    );


    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());


    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
         // Production modified due to RST correction
         G2 / nut(),
                (this->c1_/this->a1_)*this->betaStar_*omega_()
               *max(this->a1_*omega_(), this->b1_*F23()*sqrt(S2()))
            )
        // Production modified due to k-equation correction
          + alpha()*rho()*gamma*kDeficit_/nut()*(xi_ * usekDeficit_) //Removed because of instability issues
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
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
        
        gradk_ = fvc::grad(k_);

    }


    // --------------- Turbulent kinetic energy equation ----------------
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), k_)
     ==
    // Production modified due to RST correction
        alpha()*rho()*this->Pk(G2)  //this->Pk(G)
    // Production modified due to k-equation correction
      + alpha()*rho()*kDeficit_()*(xi_ * usekDeficit_)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), k_)
      + this->kSource()
      + fvOptions(alpha, rho, k_)
    );


    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    this->correctNut(S2);

    // * * * * * * * * * * * * * *  Computing (k) equation terms * * * * *** * * * * * * * * //

    // Production term for k equation
    Pk_.ref() = (alpha()*rho()*this->Pk(G2)-(2.0/3.0)*alpha()*rho()*divU*k_())/rho();

    // Destruction term of k equation
    Dk_.ref() = (alpha()*rho()*this->epsilonByk(F1, tgradU())*k_())/rho();

    tgradU.clear();

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

