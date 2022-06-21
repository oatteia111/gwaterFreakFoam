Info<<"start 2phasefields"<<endl;
Se = min((sw-sw_min)/(sw_max-sw_min),1.-SMALL);//Info<<"Se "<<Se<<endl;
volScalarField krg ("krg", pow(1-Se,0.5) * pow(1-pow(Se,1/m_vg),2*m_vg) );//Info<<"krg "<<krg<<endl;
surfaceScalarField krgf ("krgf",fvc::interpolate(krg,"krg"));
//dimensionedScalar pc0 ("pc0",dimensionSet(1,-1,-2,0,0,0,0));
dimensionedScalar pc0 = mag(g)*rhow/alpha_vg; Info<<"pc0 "<<pc0<<endl;// alphavg 1/m -> Pa
//pc0.dimensions().reset(dimless);
//volScalarField pc ("pc", unity);
//pc.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));
volScalarField dpcdS ("dpcdS", - 1/(n_vg*m_vg) * pc0/(sw_max - sw_min) * pow(pow(Se,-1/m_vg)-1,(1/n_vg)-1) * pow(Se,-(1+m_vg)/m_vg) );//Info<<"dpcdS "<<dpcdS<<endl;
volScalarField dkrgdS ("dkrgdS", unity);

	Info << "pc0 dim "<<pc0.dimensions()<<" val "<<pc0<<"dpcdS dim" << dpcdS.dimensions() << endl;
	dpcdS.dimensions().reset(dimensionSet(1,-1,-2,0,0,0,0));

//#include "createPhi.H"
	// phase
	autoPtr<incompressiblePhase> phaseg = incompressiblePhase::New(mesh,transportProperties,"g");
	volVectorField& Ug = phaseg->U();
	const dimensionedScalar& rhog = phaseg->rho();
	const dimensionedScalar& mug = phaseg->mu();
	phaseg->phi().writeOpt()=IOobject::NO_WRITE;
	dimensionedScalar Mmu(muw/mug);//viscosity ratio
	U = Ug + Uw;

Info<<"end phase g"<<endl;
//- mobility computation 
surfaceScalarField Mwf ("Mwf",Kf*krwf/muw);
surfaceScalarField Lwf ("Lwf",Kf*rhow*krwf/muw);
surfaceScalarField Fwf ("Fwf",(krwf/muw) / ( (krgf/mug) + (krwf/muw) ) );
if (activate2phase == 1) // 2 phase
{
	b.dimensions().reset(dimensionSet(-1,1,2,0,0,0,0));
	Mf.dimensions().reset(dimensionSet(-1,3,1,0,0,0,0)); // 
	Mf = Mwf + Kf*krgf/mug;Info<<"Mf 0 1 "<<Mf[0]<<" "<<Mf[1]<<endl;
	Lf.dimensions().reset(dimensionSet(0,0,1,0,0,0,0));
	Lf = Lwf + Kf*rhog*krgf/mug;//Info<<"Lf "<<Lf[10]<<endl;
}
else 
	{
		volVectorField Ug = Uw*0;
	}
Info<<"end Mf"<<endl;
//- compute fluxes
surfaceVectorField gradpc ("gradpc", fvc::interpolate(dpcdS*fvc::grad(sw),"pc") );
surfaceScalarField phiPc ("phiPc", (Mwf * gradpc) & mesh.Sf() );
phiGr = (Lf * g) & mesh.Sf();Info<<"phiGr 0 1 "<<phiGr[0]<<" "<<phiGr[1]<<endl;
surfaceScalarField phi ("phi", phiw);
Info<<"units phiGr"<<phiGr.dimensions()<<"units gradpc"<<gradpc.dimensions()<<"units phiPc"<<phiPc.dimensions()<<endl;
