if(parallelSolve=="yes")
{
	IOdictionary decomposeParDict
	(
		IOobject
		(
			"decomposeParDict",
			runTime.system(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	scalar NOS (readScalar(decomposeParDict.lookup("numberOfSubdomains")));
	
	std::stringstream call_line;
    call_line << "decomposePar -force && mpirun -np " << NOS << " rhoSimpleFoamPen -parallel";
    system(call_line.str().c_str());

    std::stringstream call_line2;
    call_line2 << "reconstructPar -latestTime";
    system(call_line2.str().c_str());
}

else system("rhoSimpleFoamPen");

system("echo Primal calculations completed!");

Info << "\nOptimization iteration " << optCont << "!" << endl;

// Read adjoint variables from last opt iteration
volScalarField psiCRead (IOobject("psiC", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volVectorField psiURead (IOobject("psiU", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volScalarField psiERead (IOobject("psiE", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);

// Read last primal solve iteration from objectiveFunction.txt
std::ifstream objFuncFile;
string It;
string JVal;

objFuncFile.open("objectiveFunction.txt");
while(!objFuncFile.eof())
{
	getline(objFuncFile, It, ' ');
	getline(objFuncFile, JVal, '\n');
}
objFuncFile.close();

int lastRunTime = std::stoi(It);
dimensionedScalar lastTime("lastTime", dimTime, lastRunTime); // Last time from primal solve
Foam::label labelDKW = 0;                                     // Label needed to use setTime (don't know why it is necessary)
runTime.setTime(lastTime, labelDKW);                          // Go to last time from primal solve

// Read primal variables from last primal solve
volScalarField pRead        (IOobject("p",         runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volVectorField URead        (IOobject("U",         runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volScalarField TRead        (IOobject("T",         runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volScalarField rhoRead      (IOobject("rho",       runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), thermo.rho());
surfaceScalarField phiRead  (IOobject("phi",       runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volScalarField rhoDesRead   (IOobject("rhoDes",    runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
volScalarField alphaPenRead (IOobject("alphaPen",  runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);

// Assign all read values
p            = pRead;
U            = URead;
thermo.he()  = thermo.Cp()*TRead; // TODO: read true enthalpy field and not Cp*T
//thermo.rho() = rhoRead; // o problema atual está aqui
phi          = phiRead;
rhoDes       = rhoDesRead;
alphaPen     = alphaPenRead;
psiC         = psiCRead;
psiU         = psiURead;
psiE         = psiERead;
