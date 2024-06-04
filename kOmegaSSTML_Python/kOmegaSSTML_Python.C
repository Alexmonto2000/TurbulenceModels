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

#include "kOmegaSSTML_Python.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTML_Python<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
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
void kOmegaSSTML_Python<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTML_Python<BasicTurbulenceModel>::kOmegaSSTML_Python
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
    sigma_
    (
        IOobject(
	    "sigma",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar
        (
            "sigma", 
            dimensionSet(0,0,0,0,0,0,0), 
            1.20813608515e-37
        ) 
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
    //========================== Python initialization ========================
    
    // Start up the Python interpreter
    Py_Initialize();
    
    // Add the run directory to the python path
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    // initialize numpy array library
    import_array1();

    // Load the file model_definition.py which should be in the case directory
    pName = PyUnicode_DecodeFSDefault("model_definition");
    pModule = PyImport_Import(pName);

    // Load the model function inside model_definition.py
    model = PyObject_GetAttrString(pModule, "model");
    
    // Initialize tuple to be send to Python
    model_args = PyTuple_New(2);


    // Set num_cells equal to the number of cells in the mesh    
    num_cells = this->mesh_.cells().size();
             
    // Initialize the (num_cells x num_scalars) sized input_vals array. It is
    // 1D on purpose so the whole array is contiguous in memory.
    input_vals = new double[num_cells*num_scalars];
}

template<class BasicTurbulenceModel>
tmp<fvVectorMatrix> kOmegaSSTML_Python<BasicTurbulenceModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_Python::divDevReff()" << endl;
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
kOmegaSSTML_Python<BasicTurbulenceModel>::devRhoReff() const
{
    Info << "In: kOmegaSSTML_Python::devRhoReff()" << endl;
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
kOmegaSSTML_Python<BasicTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_Python::divDevRhoReff(U)" << endl;
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
kOmegaSSTML_Python<BasicTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    Info << "In: kOmegaSSTML_Python::divDevRhoReff(rho, U)" << endl;
    return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U)))) //
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U) //
      + fvc::div(dev(2.*this->k_*this->bijDelta_) * useRST_ * xi_)
    );
}  
  
template<class BasicTurbulenceModel>
void kOmegaSSTML_Python<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
    Info << "In: kOmegaSSTML_Python.correct()" << endl;
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

    //========================== Python interaction during run ================

    //Load the latest p_ as a variable
    word pName_ = ("p");
    tmp<volScalarField> p_ = U.db().lookupObject<volScalarField>(pName_);     

    word lName_ = ("lorentz");
    tmp<volVectorField> lorentzforce = U.db().lookupObject<volVectorField>(lName_);     
   
    tmp<volTensorField> tgradU = fvc::grad(U);
    volSymmTensorField S(symm(tgradU()) / omega_);
    volTensorField W(skew(tgradU()) / omega_);
    volScalarField S2(2*magSqr(symm(tgradU()))); 

    //Get the gradient fields of p and k
    tmp<volVectorField> tgradp = fvc::grad(p_);
    tmp<volVectorField> tgradk = fvc::grad(k_);
    tmp<volVectorField> tcurlU = fvc::curl(U); 
    
    // Loop over each mesh cell and store relevant variables in the array which
    // is passed to Python.
    forAll(k_.internalField(), id)
    {

        // First nine elements correspond to the components of the gradU tensor
        input_vals[id*num_scalars + 0] = tgradU()[id][0];
        input_vals[id*num_scalars + 1] = tgradU()[id][1];
        input_vals[id*num_scalars + 2] = tgradU()[id][2];
        input_vals[id*num_scalars + 3] = tgradU()[id][3];
        input_vals[id*num_scalars + 4] = tgradU()[id][4];
        input_vals[id*num_scalars + 5] = tgradU()[id][5];
        input_vals[id*num_scalars + 6] = tgradU()[id][6];
        input_vals[id*num_scalars + 7] = tgradU()[id][7];
        input_vals[id*num_scalars + 8] = tgradU()[id][8];
        
        // Tenth element corresponds to k
        input_vals[id*num_scalars + 9] = k_[id];
        
        // Eleventh element corresponds to omega
        input_vals[id*num_scalars + 10] = omega_[id];

        // Twelfth-fourteenth elements correspond to components of the 
        // gradp vector.
        input_vals[id*num_scalars + 11] = tgradp()[id][0];
        input_vals[id*num_scalars + 12] = tgradp()[id][1];
        input_vals[id*num_scalars + 13] = tgradp()[id][2];

        // Fifteenth-seventeenth elements correspond to components of the
        // gradk vector.
        input_vals[id*num_scalars + 14] = tgradk()[id][0];
        input_vals[id*num_scalars + 15] = tgradk()[id][1];
        input_vals[id*num_scalars + 16] = tgradk()[id][2];
        
        // Eighteenth element corresponds to nu_t
        input_vals[id*num_scalars + 17] = nut()[id]; 
        
        // Nineteenth-twenty-first elements correspond to components of the 
        // U vector.
        input_vals[id*num_scalars + 18] = U()[id][0];
        input_vals[id*num_scalars + 19] = U()[id][1];
        input_vals[id*num_scalars + 20] = U()[id][2]; 
        
        // Twenty-second element corresponds to the wall distance
        input_vals[id*num_scalars + 21] = y_[id];   

        // Twenty-third element corresponds to the viscosity
        input_vals[id*num_scalars + 22] = this->nu()().internalField()[id];

        // Twenty-fourth-twenty-sixth elements corresponds to components
        // of curlU vector.
        input_vals[id*num_scalars + 23] = tcurlU()[id][0];
        input_vals[id*num_scalars + 24] = tcurlU()[id][1];
        input_vals[id*num_scalars + 25] = tcurlU()[id][2];

        //Lorentz force vector
        input_vals[id*num_scalars + 26] = lorentzforce()[id][0];
        input_vals[id*num_scalars + 27] = lorentzforce()[id][1];
        input_vals[id*num_scalars + 28] = lorentzforce()[id][2];
                    
    }

    // Clear temporary arrays with gradients/vorticity
    tgradp.clear();
    tgradk.clear();
    tcurlU.clear();
    
    // Get the array dimensions in a format understood by Numpy
    npy_intp dim[] = {num_cells, num_scalars};

    // Convert the input_vals array to a format understood by Numpy    
    array_2d = PyArray_SimpleNewFromData(2, dim, NPY_DOUBLE, &input_vals[0]);
    
    // Set the first element of the tuple to the array to be send to Python
    PyTuple_SetItem(model_args, 0, array_2d);
    
    // Create dictionary of booleans relevant for Python and set it as the
    // second element of the tuple to be send to Python.
    boolDict = PyDict_New();
    PyDict_SetItemString(boolDict, "modelkDeficit",
                         PyLong_FromLong(modelkDeficit_));
    PyDict_SetItemString(boolDict, "modelRST", PyLong_FromLong(modelRST_));
    PyDict_SetItemString(boolDict, "useSigma", PyLong_FromLong(useSigma_));
    PyDict_SetItemString(boolDict, "modelSigma", PyLong_FromLong(modelSigma_));
    PyTuple_SetItem(model_args, 1, boolDict);
   

    // Call the model() function in Python with the input_vals in a tuple as
    // the argument. The array returned by Python is loaded into pReturn,
    // which is initialized here.
    PyObject* pReturn = PyObject_CallObject(model, model_args);
    
    // Cast pReturn to a PyArrayObject and store it in pValue
    pValue = reinterpret_cast<PyArrayObject*>(pReturn);
       
    // Check if the returned array has the expected number of rows
    // (corresponding to number of cells). If not, throw an error.
    if (PyArray_DIMS(pValue)[0] != num_cells){
        FatalError << "Number of rows (corresponding to the number of mesh "
        << "cells) returned by Python does not correspond to the number of "
        << "mesh cells sent to Python." << nl << exit(FatalError); 
    }
    
    // Check if the returned array has the expected number of columns
    // (corresponding to first column kDeficit, next six bijDelta and last
    // column sigma). If not, throw an error.
    if (PyArray_DIMS(pValue)[1] != num_return){
        FatalError << "Number of columns of the array returned by Python does"
        << " not correspond to the expected number of columns (8). The first "
        << "column should be kDeficit, the next six columns components of "
        << "bijDelta (XX, XY, XZ, YY, ZZ) and the last column sigma." 
        << nl << exit(FatalError);
    }
    
    // If kDeficit is modeled, extract it as the first column of the returned
    // array. If it is not modeled, check whether the file was succesfully
    // read in. If not, throw an error.
    if (modelkDeficit_)
    {
        forAll(kDeficit_.internalField(), id)
        {
            kDeficit_[id] = *((double*)PyArray_GETPTR2(pValue, id, 0));
        }
    }
    else if (kDeficit_[0] == 1.20813608515e-37)
    {
        FatalError << "modelkDeficit set to false, but no kDeficit file found"
        << nl << exit(FatalError);  
    }

    // If bijDelta is modeled, extract it as the next six columns of the
    // returned array. If it is not modeled, check whether the file was
    // succesfully read in. If not, throw an error.  
    if (modelRST_)
    {
        forAll(bijDelta_.internalField(), id)
        {
            bijDelta_[id][0] = *((double*)PyArray_GETPTR2(pValue, id, 1));
            bijDelta_[id][1] = *((double*)PyArray_GETPTR2(pValue, id, 2));
            bijDelta_[id][2] = *((double*)PyArray_GETPTR2(pValue, id, 3));
            bijDelta_[id][3] = *((double*)PyArray_GETPTR2(pValue, id, 4));
            bijDelta_[id][4] = *((double*)PyArray_GETPTR2(pValue, id, 5));
            bijDelta_[id][5] = *((double*)PyArray_GETPTR2(pValue, id, 6));
      }
    }
    else if (bijDelta_[0][0] == 1.20813608515e-37)
    {
        FatalError << "modelbijDelta set to false, but no bijDelta file found"
        << nl << exit(FatalError);  
    }

    // If sigma is used, either read it in from a file or calculate it.
    // If sigma is not used, set it to 1. everywhere.
    if (useSigma_)
    {
        // If modelSigma_ is true, extract the last column of pValue as sigma. 
        // If useSigma_ is false, check if the file was succesfully read in,
        // otherwise throw an error.
        if (modelSigma_)
        {
            forAll(sigma_.internalField(), id)
            {
                sigma_[id] = *((double*)PyArray_GETPTR2(pValue, id, 7));
            }              
        }
        else if (sigma_[0] == 1.20813608515e-37)
        {
            FatalError << "modelSigma set to false, but no sigma file found"
            << nl << exit(FatalError);
        }    
    }
    else
    {
        forAll(sigma_.internalField(), id)
        {
            sigma_[id] = 1.;
        }
    }

    //Free the memory of pReturn to prevent a memory leak
    Py_DECREF(pReturn);
    
    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

/************************************************************************/

    lorentzforce.clear();

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

