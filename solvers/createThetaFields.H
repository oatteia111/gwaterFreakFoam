// base field

// get transport properties
Info << "Reading unsat properties\n" << endl;
scalar n_vg(transportProperties.lookupOrDefault<scalar>("n_vg",1.5));
dimensionedScalar alpha_vg(transportProperties.lookupOrDefault("alpha_vg",dimensionedScalar("alpha_vg",dimless/dimLength,10.0)));
scalar sw_min(transportProperties.lookupOrDefault<scalar>("sw_min",0.0001));
scalar sw_max(transportProperties.lookupOrDefault<scalar>("sw_max",0.9999));
scalar m_vg = 1 - 1/n_vg;
Info<<"step 1 "<<endl;
volScalarField thetag ("thetag",unity);
volScalarField thetaw ("thetaw",unity);
surfaceScalarField phig ("phig",phiw*0);
volScalarField krw ("krw",unity);
Info<<"step 2 "<<endl;
volScalarField dkrwdS ("dkrwdS", unity);
volScalarField Se ("Se",unity);
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));
surfaceScalarField Lf = rhow*Kf*krwf/muw;
Info<<"step 3 "<<endl;
surfaceScalarField phiGr ("phiGr",(Lf * g) & mesh.Sf() );
surfaceScalarField phiPcNull ("phiPc",0*phiGr);
label iterPicard=0;
scalar resPicard = 1000.;
Info << "end general fields" << endl;

if (activateUnsat == 1) // unsaturated
{
	//create head and calculate capillary pressure
	//hp.dimensions().reset(dimless);
	Se = pow(1.0+pow(-min(hp*alpha_vg,-1e-4),n_vg),-m_vg);
	//Se = neg(hp) *pow(1.0+pow(alpha_vg*mag(hp),n_vg),-m_vg) + pos(hp+VSMALL)*(1-SMALL);
	//hp.dimensions().reset(dimLength);
	sw = sw_min + (sw_max-sw_min)*Se;
//Info << "end Se " << endl;

	//- Update effective saturation and capillary capacity
	// HERE Ch is called stor for compatibility
	stor.dimensions().reset(dimless/dimLength);
	//Ch == alpha*m*(Smax-Smin)/(1.0-m_vg)*pow(Se,1.0/m_vg)*pow(1.0-pow(Se,1.0/m_vg),m_vg) ;
	stor = alpha_vg*m_vg*(1-sw_min)/(1-m_vg)*pow(Se,1/m_vg)*pow(1-pow(Se,1/m_vg),m_vg);
	//stor.dimensions().reset(dimless/dimLength);
//Info << "end stor " << endl;

//- relative permeability (kr)
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2); //krbmax = 1
	krwf = fvc::interpolate(krw,"krw");
	dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	dkrwdS *= 1/(1 - sw_min);
	
	Info << "ending with surf fields \n" << endl;

	//- mobility and fractional flow 
	Mf = rhow*mag(g)*Kf*krwf/muw; Info << "Mf done\n" << endl;// Mf alread exists
	Lf = rhow*Kf*krwf/muw;
	//- fluxes depending on saturation
	phiGr = (Lf * g) & mesh.Sf();

	Info << "end createThetaF \n" << endl;
}
