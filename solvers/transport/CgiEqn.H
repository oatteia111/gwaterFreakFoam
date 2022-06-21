thetag = max(eps*(sw_max - sw),1e-12);
Deffg = Dg0 * ( pow(thetag,3.333) + 1e-4*pow(eps*sw,3.333) ) /pow(eps,2.)*tensor::I;
Info<<"sw "<<sw[2]<<" thetag "<<thetag[2]<<" Deffg 2 "<<Deffg[2]<<endl;

for (i=0; i<ph_gcomp;i++)
    {
		//Info << "in CgiEqn "<< i << endl;
		volScalarField& Cgi = Cg[i]();
		word nameCi ("Cg" + name(i));
		
		fvScalarMatrix CgiEqn
		(		  
            thetag * fvm::ddt(Cgi) +  Cgi * fvc::ddt(thetag)
			+ fvm::div(phig, Cgi, "div(phig,Cgi)") // used only in multiphase
            - fvm::laplacian(Deffg, Cgi, "laplacian(Deffg,Cgi)")
			==
			fvOptions(sw,Cgi)
		);

		fvOptions.constrain(CgiEqn);
		CgiEqn.solve();
    }
