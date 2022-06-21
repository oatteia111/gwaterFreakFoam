{
// diffusion+dispersion
dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
//Deffw = (alphaT * mag(vp))*tensor::I + (alphaL-alphaT) * (vp*vp)/(mag(vp)+USMALL); // Ut rpelaced by vp
Deffw = (alphaT * mag(Uw))*tensor::I + (alphaL-alphaT) * (Uw*Uw)/(mag(Uw)+USMALL);
forAll(Deffw.boundaryField(),patchi)
	{ //- Dispersion equal to zero on patch where flux is fixed
	if (Cw.boundaryField().types()[patchi] == "fixedFlux") Deffw.boundaryFieldRef()[patchi] = tensor(zero());
	}
thetaw = eps*sw;

fvScalarMatrix CwEqn
	(
		thetaw * fvm::ddt(Cw) + Cw * fvc::ddt(thetaw)
		+ fvm::div(phiw, Cw, "div(phiw,Cw)") 
		- fvm::laplacian(Deffw, Cw, "laplacian(Deffw,Cw)")
		==
		fvOptions(sw,Cw) 
	);
fvOptions.constrain(CwEqn);
CwEqn.solve();
//Info << "Ceqn solve "<<endl;
}
