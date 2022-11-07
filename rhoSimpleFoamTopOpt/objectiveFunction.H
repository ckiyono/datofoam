scalar J = 0;

if(objFunc == "entropyVar") // entropy variation (for perfect gas)
{
	dimensionedScalar pRef("pRef", dimMass*dimLength/(dimTime*dimTime*dimLength*dimLength), readScalar(thermophysicalProperties.lookup("pRef")));
	dimensionedScalar TRef("TRef", dimTemperature, readScalar(thermophysicalProperties.lookup("TRef")));
	dimensionedScalar molWeight("molWeight", dimensionSet(1,0,0,0,-1,0,0), readScalar(thermophysicalProperties.subDict("mixture").subDict("specie").lookup("molWeight")));
	dimensionedScalar RUniversal("RUniversal", dimensionSet(1,2,-2,-1,-1,0,0), (Foam::constant::thermodynamic::RR));

	dimensionedScalar Rgas = RUniversal/molWeight;
	
	sEnt = sRefField + thermo.Cp()*log(thermo.T()/TRef)-Rgas*log(p/pRef);
		
	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].type() == "patch")
			{
				J = J +
				(sEnt.boundaryField()[patchI][faceI]*rho.boundaryField()[patchI][faceI])*
				(U.boundaryField()[patchI][faceI] & mesh.Sf().boundaryField()[patchI][faceI]);
			}
		}
	}
	
	reduce (J, sumOp<scalar>());
}

if(objFunc == "impellerEff") // Isentropic efficiency
{	
	scalar isoWork  = 0;
	scalar realWork = 0;
	
	volVectorField Uabs  = U + (omegaVector^mesh.C());
	volScalarField hTot  = thermo.he() + 0.5*(Uabs & Uabs);
	
	scalar areaOut     = 0;
	scalar mOut        = 0;
	scalar pOutCalc    = 0;
	scalar hTotOutCalc = 0;
	
	scalar areaIn     = 0;
	scalar mIn        = 0;
	scalar pInCalc    = 0;
	scalar hTotInCalc = 0;
	
	scalar TInCalc = 0;
	
	faceSet pRatioOut
	(
		IOobject
		(
			"pRatioOut", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);
	
	faceSet pRatioIn
	(
		IOobject
		(
			"pRatioIn", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);

	const labelList& pRatioOutList = pRatioOut.toc();
	const labelList& pRatioInList  = pRatioIn.toc();

	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].name() == "outlet")
			{
				scalar dm = rho.boundaryField()[patchI][faceI]*(U.boundaryField()[patchI][faceI] & mesh.Sf().boundaryField()[patchI][faceI]);
				areaOut   = areaOut + mag(mesh.Sf().boundaryField()[patchI][faceI]);
				
				mOut        = mOut + dm;
				pOutCalc    = pOutCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				hTotOutCalc = hTotOutCalc + hTot.boundaryField()[patchI][faceI]*dm;
			}
			
			if(mesh.boundary()[patchI].name() == "inlet")
			{
				scalar dm = rho.boundaryField()[patchI][faceI]*(U.boundaryField()[patchI][faceI] & mesh.Sf().boundaryField()[patchI][faceI]);
				areaIn    = areaIn + mag(mesh.Sf().boundaryField()[patchI][faceI]);
				
				mIn        = mIn + mag(dm);
				pInCalc    = pInCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				hTotInCalc = hTotInCalc + hTot.boundaryField()[patchI][faceI]*mag(dm);
				
				TInCalc = TInCalc + thermo.T().boundaryField()[patchI][faceI]*mag(dm);
			}
		}
	}
	
	reduce (areaOut, sumOp<scalar>());
	reduce (mOut, sumOp<scalar>());
	reduce (pOutCalc, sumOp<scalar>());
	reduce (hTotOutCalc, sumOp<scalar>());
	
	reduce (areaIn, sumOp<scalar>());
	reduce (mIn, sumOp<scalar>());
	reduce (pInCalc, sumOp<scalar>());
	reduce (hTotInCalc, sumOp<scalar>());
	
	reduce (TInCalc, sumOp<scalar>());
	
	scalar pOut    = pOutCalc/areaOut;
	scalar hTotOut = hTotOutCalc/mOut;
	
	scalar pIn     = pInCalc/areaIn;
	scalar hTotIn  = hTotInCalc/mIn;
	
	scalar TIn  = TInCalc/mIn;
	
	scalar TOutIso = TIn*Foam::pow((pOut/pIn), (0.4/1.4));
	
	scalar cpCalc(readScalar(thermophysicalProperties.subDict("mixture").subDict("thermodynamics").lookup("Cp")));
	
	isoWork  = cpCalc*TOutIso - hTotIn;
	realWork = hTotOut - hTotIn;
	
	J = -(isoWork/realWork);
}

if(objFunc == "impellerTorque")
{	
	vector torque   = vector (0, 0, 0);
	volVectorField Uabs = U + (omegaVector^mesh.C());
	
	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].type() == "patch")
			{
				int faceID = mesh.boundary()[patchI].start() + faceI;
				
				torque = torque +
				rho.boundaryField()[patchI][faceI]*
				(mesh.Cf()[faceID]^(Uabs.boundaryField()[patchI][faceI]))*
				(Uabs.boundaryField()[patchI][faceI] & mesh.Sf().boundaryField()[patchI][faceI]);					
			}
		}
	}

	reduce (torque, sumOp<vector>());	

	J = torque.component(2);
}
