/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::mixingPrecipitation

Description
    Precipitation model based on the Yang et al. (2021) surface rate law

        rDot = kPrec*(Omega - 1)^nPrec,  for Omega > 1

    This solver reads separate dissolved cation and anion concentration fields
    and computes a concentration-based saturation ratio:

        Omega = Ccation*Canion/Ksp

    The surface molar rate rDot is converted to a volumetric mass source with
    the same area-density and interface-mask treatment used by the existing
    precipitate model.

Usage

    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to solid)
        {
            type                mixingPrecipitation;
            kPrec               1e-6;
            nPrec               1.0;
            Ksp                 4.0;
            Mv                  1000.0;
            alphaSolidMin       0.99;
            smoothSurface       true;
            smoothAreaDensity   10;
            reactiveWallPatches (walls frontAndBack);
            writeSaturationFields true;
        }
    );
    \endverbatim

    where:
        kPrec is in mol/m2/s
        nPrec is dimensionless
        Ksp is in (mol/m3)^2 for the concentration-based approximation
        Mv follows the existing case convention and is read in g/mol

SourceFiles
    mixingPrecipitation.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_mixingPrecipitation_H
#define meltingEvaporationModels_mixingPrecipitation_H

#include "InterfaceCompositionModel.H"

namespace Foam
{

class phasePair;

namespace meltingEvaporationModels
{

template<class Thermo, class OtherThermo>
class mixingPrecipitation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Surface molar precipitation coefficient from Eq. (7)
        dimensionedScalar kPrec_;

        //- Kept for interface compatibility with the existing base class
        const dimensionedScalar Tactivate_;

        //- Precipitate molar mass, using the existing g/mol case convention
        dimensionedScalar Mv_;

        //- Reaction order in (Omega - 1)^n
        scalar nPrec_;

        //- Solubility product for concentration-based Omega
        scalar Ksp_;

        //- Concentration field names used to compute Omega
        word cationFieldName_;
        word anionFieldName_;

        //- 'To' phase maximum value for the mass transfer
        scalar alphaMax_;

        //- 'To' phase minimum value for the mass transfer
        scalar alphaMin_;

        //- 'To' phase minimum value for the surface identification
        scalar alphaSolidMin_;

        //- Alpha maximum for the rest of phases
        scalar alphaRestMax_;

        //- Surface area-density option
        bool smoothSurface_;

        //- User-input smooth surface area-density constant [1/m]
        scalar smoothAreaDensity_;

        //- Boundary patches treated as pre-existing solid nucleation surfaces
        wordList reactiveWallPatches_;

        //- Optional write of Omega and SI fields on output times
        bool writeSaturationFields_;


public:

    TypeName("mixingPrecipitation");


    mixingPrecipitation
    (
        const dictionary& dict,
        const phasePair& pair
    );


    virtual ~mixingPrecipitation() = default;


    virtual tmp<volScalarField> Kexp
    (
        label variable,
        const volScalarField& field
    );

    virtual const dimensionedScalar& Tactivate() const;
};

} // End namespace meltingEvaporationModels
} // End namespace Foam

#ifdef NoRepository
#   include "mixingPrecipitation.C"
#endif

#endif

// ************************************************************************* //
