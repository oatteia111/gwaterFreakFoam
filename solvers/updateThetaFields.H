//- update saturation 
	//hp = p/rho/mag(g);
	//hp.dimensions().reset(dimless);
	Se =  pow(1.0+pow(-min(hp*alpha_vg,-1e-4),n_vg),-m_vg);
	//Se =  neg(hp) *pow(1.0+pow(alpha_vg*mag(hp),n_vg),-m_vg) + pos(hp+VSMALL)*(1-SMALL);
	Info << "n "<<n_vg<<" alpha "<<alpha_vg<<" m_vg "<<m_vg<<" Se 10 "<<Se[10]<<" Se 20 "<<Se[20]<<" Se 30 "<<Se[30]<<endl;
	//hp.dimensions().reset(dimLength);
	sw = sw_min + (sw_max-sw_min)*Se;
	//forAll(theta,i){Info<<"theta "<<theta[i]<<endl;}

	//- Update effective saturation and capillary capacity (stor=Ch)
	//stor.dimensions().reset(dimless/dimLength);
	stor.dimensions().reset(dimless/dimLength);
	stor = eps*alpha_vg*m_vg*(1.-sw_min)/(1.-m_vg)*pow(Se,1./m_vg)*pow(1.-pow(Se,1./m_vg),m_vg) ;
	//Info<<"stor "<<stor[10]<<endl;

//- relative permeability computation
	//Se==(Sb-Smin)/(Smax-Smin);
	krw = 1 * pow(Se,0.5) * pow(1-pow(1-pow(Se,1/m_vg),m_vg),2); //krmax = 1
	//Info<<"kr "<<kr[10]<<endl;

	dkrwdS = 0.5 * (1-pow((1-pow(Se,1/m_vg)),m_vg));
	dkrwdS *= ( 4 * pow(Se,1/m_vg-1/2) * pow( (1-pow(Se,1/m_vg)) , m_vg-1)) - ( pow((1-pow(Se,1/m_vg)),m_vg) -1) / pow(Se,0.5);
	dkrwdS *= 1/(1 - sw_min);
	//Info<<"dkr "<<dkrdS[10]<<endl;

krwf = fvc::interpolate(krw,"krw");

//- mobility computation 
Mf = krwf*Kf*rhow*mag(g)/muw;
	//Info<<"Mf "<<Mf[10]<<endl;
Lf = rhow*Kf*krwf/muw;
	//Info<<"Lf "<<Lf[10]<<endl;


//- compute fluxes
phiGr = (Lf * g) & mesh.Sf();
	//Info<<"PhiG "<<phiG[10]<<endl;
