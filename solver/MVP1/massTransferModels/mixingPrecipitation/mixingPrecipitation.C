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

\*---------------------------------------------------------------------------*/

#include "mixingPrecipitation.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::mixingPrecipitation<Thermo, OtherThermo>
::mixingPrecipitation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    kPrec_("kPrec", dimMoles/dimArea/dimTime, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    Mv_("Mv", dimMass/dimMoles, -1, dict),
    nPrec_(dict.lookupOrDefault<scalar>("nPrec", 1.0)),
    Ksp_
    (
        dict.lookupOrDefault<scalar>
        (
            "Ksp",
            dict.lookupOrDefault<scalar>("Kbarite", -1.0)
        )
    ),
    cationFieldName_(dict.lookupOrDefault<word>("cationField", "Ccation")),
    anionFieldName_(dict.lookupOrDefault<word>("anionField", "Canion")),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.5)),
    alphaSolidMin_(dict.lookupOrDefault<scalar>("alphaSolidMin", 0.5)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01)),
    smoothSurface_(dict.lookupOrDefault<bool>("smoothSurface", false)),
    smoothAreaDensity_(dict.lookupOrDefault<scalar>("smoothAreaDensity", 1.0)),
    reactiveWallPatches_
    (
        dict.lookupOrDefault<wordList>("reactiveWallPatches", wordList())
    ),
    writeSaturationFields_
    (
        dict.lookupOrDefault<bool>("writeSaturationFields", false)
    )
{
    if (this->transferSpecie() != "none")
    {
        word fullSpeciesName = this->transferSpecie();
        const auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));

        const typename OtherThermo::thermoType& toThermo =
            this->getLocalThermo(speciesName, this->toThermo_);

        // Keep the existing g/mol convention used by the case dictionaries.
        Mv_.value() = toThermo.W();
    }

    if (Mv_.value() == -1)
    {
        FatalErrorInFunction
            << "Please provide the precipitate molar weight Mv [g/mol]"
            << abort(FatalError);
    }

    if (kPrec_.value() < 0)
    {
        FatalErrorInFunction
            << "kPrec must be non-negative"
            << abort(FatalError);
    }

    if (nPrec_ < 0)
    {
        FatalErrorInFunction
            << "nPrec must be non-negative"
            << abort(FatalError);
    }

    if (Ksp_ <= 0)
    {
        FatalErrorInFunction
            << "Ksp must be positive"
            << abort(FatalError);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::mixingPrecipitation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field)
{
    const volScalarField& to = this->pair().to();
    const volScalarField& from = this->pair().from();
    const fvMesh& mesh = this->mesh_;
    const volScalarField& Ccation =
        mesh.lookupObject<volScalarField>(cationFieldName_).oldTime();
    const volScalarField& Canion =
        mesh.lookupObject<volScalarField>(anionFieldName_).oldTime();

    const dimensionedScalar molarMassKgPerMol
    (
        "molarMassKgPerMol",
        dimMass/dimMoles,
        Mv_.value()*1e-3
    );

    const volVectorField gradFrom(fvc::grad(from));
    const volVectorField gradTo(fvc::grad(to));

    const volScalarField areaDensitySmooth
    (
        IOobject
        (
            "areaDensitySmooth",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "areaDensitySmoothInput",
            dimless/dimLength,
            smoothAreaDensity_
        )
    );

    const volScalarField areaDensityGrad
    (
        "areaDensityGrad",
        2*mag(gradFrom)
    );

    const volScalarField gradAlphaf(gradFrom & gradTo);

    volScalarField Cmask("Cmask", from*0.0);
    Cmask.correctBoundaryConditions();

    forAll(Cmask, celli)
    {
        if (gradAlphaf[celli] < 0)
        {
            if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
            {
                const scalar alphaRes = 1.0 - from[celli] - to[celli];

                if (alphaRes < alphaRestMax_)
                {
                    Cmask[celli] = 1.0;
                }
            }
        }

        bool flag = false;

        forAll(mesh.cellCells()[celli], cellj)
        {
            if (to[mesh.cellCells()[celli][cellj]] > alphaSolidMin_)
            {
                flag = true;
            }
        }

        Cmask[celli] = flag ? 1.0 : 0.0;
    }

    forAll(mesh.boundaryMesh(), patchI)
    {
        const fvPatchScalarField& pf = to.boundaryField()[patchI];
        const labelList& faceCells = pf.patch().faceCells();
        bool reactiveWallPatch = false;

        forAll(reactiveWallPatches_, reactivePatchI)
        {
            if (pf.patch().name() == reactiveWallPatches_[reactivePatchI])
            {
                reactiveWallPatch = true;
                break;
            }
        }

        if (reactiveWallPatch)
        {
            forAll(faceCells, faceI)
            {
                Cmask[faceCells[faceI]] = 1.0;
            }
        }

        if (pf.coupled())
        {
            const scalarField neighbors = pf.patchNeighbourField();

            forAll(faceCells, faceI)
            {
                if (neighbors[faceI] > alphaSolidMin_)
                {
                    Cmask[faceCells[faceI]] = 1.0;
                }
            }
        }
    }

    volScalarField Omega
    (
        IOobject
        (
            "Omega",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    volScalarField SI
    (
        IOobject
        (
            "SI",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    volScalarField surfaceRateMol
    (
        IOobject
        (
            "surfaceRateMol",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMoles/dimArea/dimTime, Zero)
    );

    forAll(mesh.C(), cellI)
    {
        const scalar cationCell = max(Ccation[cellI], 0.0);
        const scalar anionCell = max(Canion[cellI], 0.0);
        const scalar omega =
            max(cationCell*anionCell/(Ksp_ + VSMALL), 0.0);

        Omega[cellI] = omega;
        SI[cellI] = Foam::log10(max(omega, SMALL));
        surfaceRateMol[cellI] =
            kPrec_.value()*Foam::pow(max(omega - 1.0, 0.0), nPrec_);
    }

    volScalarField massFluxPrec
    (
        "massFluxPrec",
        surfaceRateMol*molarMassKgPerMol*Cmask*pos(from - 0.01)
    );

    if (mesh.time().outputTime())
    {
        areaDensityGrad.write();
        Cmask.write();

        if (writeSaturationFields_)
        {
            Omega.write();
            SI.write();
        }
    }

    const scalar dt = mesh.time().deltaTValue();

    if (smoothSurface_)
    {
        forAll(mesh.C(), cellI)
        {
            const scalar availableMass =
                min(max(Ccation[cellI], 0.0), max(Canion[cellI], 0.0))
               *molarMassKgPerMol.value();

            if
            (
                massFluxPrec[cellI]*areaDensitySmooth[cellI]*dt
              > availableMass
            )
            {
                massFluxPrec[cellI] =
                    availableMass
                   /(dt*areaDensitySmooth[cellI] + VSMALL);
            }
        }
    }
    else
    {
        forAll(mesh.C(), cellI)
        {
            const scalar availableMass =
                min(max(Ccation[cellI], 0.0), max(Canion[cellI], 0.0))
               *molarMassKgPerMol.value();

            if
            (
                massFluxPrec[cellI]*areaDensityGrad[cellI]*dt
              > availableMass
            )
            {
                massFluxPrec[cellI] =
                    availableMass
                   /(dt*areaDensityGrad[cellI] + VSMALL);
            }
        }
    }

    if (smoothSurface_)
    {
        return massFluxPrec*areaDensitySmooth;
    }

    return massFluxPrec*areaDensityGrad;
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::mixingPrecipitation<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
