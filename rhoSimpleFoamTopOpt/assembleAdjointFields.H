forAll(mesh.C(), row)
{
	psiU[row][0] = adjSolve[3*row];
	psiU[row][1] = adjSolve[3*row + 1];
	psiU[row][2] = adjSolve[3*row + 2];;
	psiC[row]    = adjSolve[3*nCells + row];
	psiE[row]    = adjSolve[(3 + 1)*nCells + row];
}

adjSolve = -dRdvar*adjSolve;

// Note: in the sensitivities assignement, the last term
//       inside parentheses is relative to the Riesz
//       Representation Theorem. For Finite Volumes, the correct
//       would be only divide by the cell volume. The
//       multiplication by the average cell volume is only to
//       not change scale so much.

forAll(mesh.C(), cellI) Sensitivities[cellI] = adjSolve[cellI]*((gSum(mesh.V())/nCells)/mesh.V()[cellI]);
// forAll(mesh.C(), cellI) Sensitivities[cellI] = adjSolve[cellI];

psiC.write();
psiU.write();
psiE.write();
Sensitivities.write();

Info << "Sensitivities calculation completed!\n" << endl;

if((optConst == "velocProfile") or (optConst == "pRatio"))
{
	volScalarField constSens
	(
		IOobject
		(
			"constSens",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0),0)
	);
	adjConstSolve = -dRdvar*adjConstSolve;

	forAll(mesh.C(), cellI) constSens[cellI] = adjConstSolve[cellI]*((gSum(mesh.V())/nCells)/mesh.V()[cellI]);
	// forAll(mesh.C(), cellI) constSens[cellI] = adjConstSolve[cellI];
	constSens.write();
}
