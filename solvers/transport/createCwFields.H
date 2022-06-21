dimensionedScalar alphaL(transportProperties.lookupOrDefault("alphaL",dimensionedScalar("alphaL",dimLength,0.0)));
dimensionedScalar alphaT(transportProperties.lookupOrDefault("alphaT",dimensionedScalar("alphaT",dimLength,0.0)));
dimensionedScalar lbdaw(transportProperties.lookupOrDefault("lbdaw",dimensionedScalar("lbdaw",dimless/dimTime,0.0)));

Info<< "Reading field Cw\n" << endl;
wordList bcTypes (Kh.boundaryField().size(),zeroGradientFvPatchScalarField::typeName);

volScalarField Cw
(
    IOobject("Cw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE),mesh,dimensionedScalar("",dimMass/dimVol,0),bcTypes
);

volTensorField Deffw
(
    IOobject("Deffw",runTime.timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::NO_WRITE),mesh,dimensionedTensor("",dimArea/dimTime,Zero) 
);

Info<< "All C read\n" << endl;
