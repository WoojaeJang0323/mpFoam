    Info<< "Reading personlized input parameters" << endl;
    IOdictionary inputParameters
(
    IOobject
    (
        "inputParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

	    // Molar weight of precipitated species
	    dimensionedScalar molW
	    (
		    "molW",
		    dimMass/dimMoles,
		    inputParameters
	    );

        const Switch useTwoSpecies
        (
            inputParameters.lookupOrDefault<Switch>("useTwoSpecies", false)
        );

	    // Diffusion coefficient of acid in fluid phase
	    dimensionedScalar diffInLiquid
	    (
	    	"diffInLiquid",
		    dimLength*dimLength/dimTime,
		    inputParameters
	    );

        dimensionedScalar diffCationInLiquid
        (
            "diffCationInLiquid",
            dimLength*dimLength/dimTime,
            inputParameters.lookupOrDefault<scalar>
            (
                "diffCationInLiquid",
                diffInLiquid.value()
            )
        );

        dimensionedScalar diffAnionInLiquid
        (
            "diffAnionInLiquid",
            dimLength*dimLength/dimTime,
            inputParameters.lookupOrDefault<scalar>
            (
                "diffAnionInLiquid",
                diffInLiquid.value()
            )
        );

        dimensionedScalar nuCation
        (
            "nuCation",
            dimless,
            inputParameters.lookupOrDefault<scalar>("nuCation", 1.0)
        );

        dimensionedScalar nuAnion
        (
            "nuAnion",
            dimless,
            inputParameters.lookupOrDefault<scalar>("nuAnion", 1.0)
        );

    Info<< "The molar weight of solid is set as: \n" << molW << endl;

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Note: construct T to be around before the thermos. The thermos will
    //       not update T.
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	    Info<< "Reading field C\n" << endl;
	    volScalarField C 
	    (
	        IOobject
	        (
	            "C",
	            runTime.timeName(),
	            mesh,
	            IOobject::READ_IF_PRESENT,
	            IOobject::AUTO_WRITE
	        ),
	        mesh,
            dimensionedScalar("CDefault", dimMoles/dimVolume, Zero)
	    );

        Info<< "Reading field Ccation\n" << endl;
        volScalarField Ccation
        (
            IOobject
            (
                "Ccation",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("CcationDefault", dimMoles/dimVolume, Zero)
        );

        Info<< "Reading field Canion\n" << endl;
        volScalarField Canion
        (
            IOobject
            (
                "Canion",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("CanionDefault", dimMoles/dimVolume, Zero)
        );
	
	    volScalarField SuOut 
	    (
	        IOobject
        (
            "SuOut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	        ),
	        mesh
	    );

        volScalarField SuCationOut
        (
            IOobject
            (
                "SuCationOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("SuCationOutDefault", dimMoles/dimVolume/dimTime, Zero)
        );

        volScalarField SuAnionOut
        (
            IOobject
            (
                "SuAnionOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("SuAnionOutDefault", dimMoles/dimVolume/dimTime, Zero)
        );

	    Info<< "Reading field diffusion coefficient for C field\n" << endl;
	    volScalarField diffCo
    (
        IOobject
	        (
	            "diffCo",
	            runTime.timeName(),
	            mesh,
	            IOobject::READ_IF_PRESENT,
	            IOobject::AUTO_WRITE
	        ),
	        mesh,
            dimensionedScalar("diffCoDefault", dimLength*dimLength/dimTime, Zero)
	    );

        volScalarField diffCoCation
        (
            IOobject
            (
                "diffCoCation",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "diffCoCationDefault",
                dimLength*dimLength/dimTime,
                Zero
            )
        );

        volScalarField diffCoAnion
        (
            IOobject
            (
                "diffCoAnion",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "diffCoAnionDefault",
                dimLength*dimLength/dimTime,
                Zero
            )
        );

    volScalarField areaDensityGrad 
    (
         IOobject
         (
             "areaDensityGrad",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         mesh
    );

    Info<< "Calculating field g.h\n" << endl;
    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"
    #include "gh.H"

    volScalarField  p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );


    Info<< "Creating multiphaseSystem\n" << endl;
    autoPtr<multiphaseSystem> fluidPtr = multiphaseSystem::New(mesh);

    multiphaseSystem& fluid = fluidPtr();

    if (!fluid.incompressible())
    {
        FatalError << "One or more phases are not incompressible. " << nl
            << "This is a incompressible solver." << abort(FatalError);
    }


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.rho()
    );
    rho.oldTime();

    // Update p using fluid.rho()
    p = p_rgh + rho*gh;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    // Mass flux
    surfaceScalarField& rhoPhi = fluid.rhoPhi();
    

    // Construct incompressible turbulence model
    autoPtr<CompressibleTurbulenceModel<multiphaseSystem>> turbulence
    (
        CompressibleTurbulenceModel<multiphaseSystem>::New
        (
            rho,
            U,
            rhoPhi,
            fluid
        )
    );

    // Creating radiation model
//    autoPtr<radiation::radiationModel> radiation
//    (
//        radiation::radiationModel::New(T)
//    );

    Info<< "Calculating field rhoCp\n" << endl;
    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.rho()*fluid.Cp()
    );
    rhoCp.oldTime();


    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
