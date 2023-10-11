/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "tractionDisplacementFvPatchVectorField.H"
#include "solidDisplacementThermo.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::tractionDisplacementFvPatchVectorField::updateCoeffs
(
    const Type& pressure
)
{
    const dictionary& mechanicalProperties =
        db().lookupObject<IOdictionary>("mechanicalProperties");
       
        
    const fvPatchField<scalar>& xh =
        patch().lookupPatchField<volScalarField, scalar>("xh");
    
    scalar rhoE(readScalar(mechanicalProperties.lookup("rhoE")));
    scalar Po(readScalar(mechanicalProperties.lookup("Po")));
    scalar rho(readScalar(mechanicalProperties.lookup("rho")));
    scalar E(rhoE/rho);
    scalar Emin(E*1e-4);
    Switch planeStress(mechanicalProperties.lookup("planeStress"));

    scalarField mu((Emin+xh*(E-Emin))/(2.0*(1.0 + Po)));
    scalarField lambda((Emin+xh*(E-Emin))*Po/((1.0 + Po)*(1.0 - 2.0*Po)));
    if (planeStress)
    {
        lambda=(Emin+xh*(E-Emin))*Po/((1.0 + Po)*(1.0 - Po));
    }
    scalarField twoMuLambda(2*mu + lambda);

    vectorField n(patch().nf());
    //vectorField mm(patch().Cf());
    //vector m(0,0,1);
    const fvPatchField<symmTensor>& sigmaD =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

    gradient() =
    (
        //(traction_ - pressure_*(n^mm))/rho
        (traction_ - pressure*n)/rho
      + twoMuLambda*fvPatchField<vector>::snGrad() - (n & sigmaD)
    )/twoMuLambda;
    fixedGradientFvPatchVectorField::updateCoeffs();
}


// ************************************************************************* //
