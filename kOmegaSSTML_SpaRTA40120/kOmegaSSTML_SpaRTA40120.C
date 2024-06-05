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

#include "kOmegaSSTML_SpaRTA40120.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
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
void kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::kOmegaSSTML_SpaRTA40120
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
tmp<fvVectorMatrix> kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA40120::divDevReff()" << endl;
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
kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::devRhoReff() const
{
    Info << "In: kOmegaSSTML_SpaRTA40120::devRhoReff()" << endl;
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
kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA40120::divDevRhoReff(U)" << endl;
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
kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_SpaRTA40120::divDevRhoReff(rho, U)" << endl;
    return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U)))) //
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U) //
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)
    );
}  
  
template<class BasicTurbulenceModel>
void kOmegaSSTML_SpaRTA40120<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    Info << "In: kOmegaSSTML_SpaRTA40120.correct()" << endl;
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
    word lName_ = ("lorentz");
    tmp<volVectorField> lorentz = U.db().lookupObject<volVectorField>(lName_);    
        
    tmp<volTensorField> tgradU = fvc::grad(U);
    volSymmTensorField S(symm(tgradU()) / omega_);
    volTensorField W(skew(tgradU()) / omega_);
    volScalarField S2(2*magSqr(symm(tgradU()))); 
    volScalarField tau_s( 1 / sqrt(magSqr(tgradU())));
    volSymmTensorField S_s(symm(tgradU()) * tau_s);
    volTensorField W_s(skew(tgradU()) * tau_s);
    /*******************************Feature Addition**************************/
    //add features for models
    // volSymmTensorField T1(S);
    // volSymmTensorField T2((symm((S & W) - (W & S))));
    // volSymmTensorField T3((symm((S & S) - 0.3333*tr(S & S)*I)));
    // volSymmTensorField T4((symm((W & W) - 0.3333*tr(W & W)*I)));
    volSymmTensorField T6((symm(W_s & W_s & S_s + S_s & W_s & W_s - 0.6666667* tr(S_s & W_s & W_s)*I)));

    gradk_ = fvc::grad(k_);

    // Basis vectors
    vector ex(1, 0, 0);
    vector ey(0, 1, 0);
    vector ez(0, 0, 1);

    // Compute the cross products with each basis vector using OpenFOAM's built-in functions
    tmp<volVectorField> crossEx = gradk_ ^ ex;
    tmp<volVectorField> crossEy = gradk_ ^ ey;
    tmp<volVectorField> crossEz = gradk_ ^ ez;

    // Initialize the result tensor field with an IOobject
    tmp<volTensorField> I_cross_gradk = W_s * dimensionedScalar("small", dimensionSet(0,1,-2,0,0,0,0), 1) ;
    // Stack the cross product vector fields into the tensor field
    // I_cross_gradk().component(0) = crossEx().component(0);
    // I_cross_gradk().component(1) = crossEy().component(0);
    // I_cross_gradk().component(2) = crossEz().component(0);
    for (int i = 0; i < 3; ++i) {
        I_cross_gradk().component(0+i*3) = crossEx().component(i);
        I_cross_gradk().component(1+i*3) = crossEy().component(i);
        I_cross_gradk().component(2+i*3) = crossEz().component(i);
    }
    crossEx.clear();
    crossEy.clear();
    crossEz.clear();

    tmp<volVectorField> crossExg = lorentz() ^ ex;
    tmp<volVectorField> crossEyg = lorentz() ^ ey;
    tmp<volVectorField> crossEzg = lorentz() ^ ez;

    // Initialize the result tensor field
    tmp<volTensorField> I_cross_lorentz = W_s * dimensionedScalar("small", dimensionSet(1,-2,-2,0,0,0,0), 1)  ;

    // Stack the cross product vector fields into the tensor field
    // I_cross_lorentz().component(0) = crossExg().component(0);
    // I_cross_lorentz().component(1) = crossEyg().component(0);
    // I_cross_lorentz().component(2) = crossEzg().component(0);
    for (int i = 0; i < 3; ++i) {
        I_cross_lorentz().component(0+i*3) = crossExg().component(i);
        I_cross_lorentz().component(1+i*3) = crossEyg().component(i);
        I_cross_lorentz().component(2+i*3) = crossEzg().component(i);
    }

    crossExg.clear();
    crossEyg.clear();
    crossEzg.clear();

    volTensorField Ak(-(tau_s /sqrt(k_)) * I_cross_gradk());
    volTensorField Alf(-(tau_s /(rhol * sqrt(k_))) * I_cross_lorentz());
    
    I_cross_gradk.clear();
    I_cross_lorentz.clear();
    // I_cross_gradk.clear();
    // I_cross_lorentz.clear();

    epsilon_ = (omega_*k_);
    volScalarField G1_s(2 * k_ * S_s && tgradU());
    volScalarField G6_s(2 * k_ * T6 && tgradU());
    
    volScalarField q_Re( 2 - min(sqrt(k_)*y_ / (50 * nu),dimensionedScalar("small", dimensionSet(0,0,0,0,0,0,0), 2)));
    volScalarField Alf2(tr(Alf & Alf));
    volScalarField Ak2(tr(Ak & Ak));
    volScalarField AlfAkS_s(tr(Alf & Ak & S_s));
    volScalarField Ak2WS_s(tr(Ak & Ak & W_s & S_s));
/************************************************************************/


/*******************************Model Addition**************************/
    //Model kDeficit_
    volScalarField classifier(0.5 + 0.5*tanh(10000*(q_Re - 0.75)));
    kDeficit_ = max(classifier * ((-1.6724144 * pow(Ak2/0.66340377 , 2) * (G6_s/0.016290889) +
     9.7697019 * (AlfAkS_s/0.4357038)*(G1_s/0.032581755) + 0.0068832524 * exp((q_Re/0.41004645))*(G1_s/0.032581755) +
      -120.9269 * ((Alf2/31.764739)/ (1 + pow(Alf2/31.764739, 2)))*(G1_s/0.032581755) + 0.90196371 * ((q_Re/0.41004645)/ (1 + pow(q_Re/0.41004645 , 2)))*(G1_s/0.032581755) +
       -3.0523038 * sqrt(pow(AlfAkS_s/0.4357038 , 2))*(G1_s/0.032581755)) * 0.0060813918),
     dimensionedScalar("small", dimensionSet(0,2,-3,0,0,0,0), 0));
/************************************************************************/

    lorentz.clear();

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
          // + alpha()*rho()*gamma*kDeficit_/nut()*(xi_ * usekDeficit_) //Removed because of instability issues
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

