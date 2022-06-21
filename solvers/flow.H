Uw.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
	if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
		{ phiw.boundaryFieldRef()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
} 

if (activateUnsat == 1) // unsaturated
{
   Info <<"in unsat "<<endl;
	iterPicard = 0; resPicard = 1000.;
	while (resPicard > tolPicard)
	{
		iterPicard++;
		#include "updateThetaFields.H"
		Info <<" iter picard "<<iterPicard<<" sw 90 "<<sw[90]<<" krf 90 "<<krwf[90]<<" Mf 90 "<<Mf[90]<<endl;
		#include "hpEqn.H"
		if (iterPicard == maxIterPicard)
		{
			Warning() <<  " Max iteration reached in Picard loop" << endl;
			break;
		}
	}
	volScalarField dsw_tmp = mag(sw-sw.oldTime());
	dsw = gMax(dsw_tmp);
	dswdTmax = dsw/runTime.deltaTValue();
}
else if (activate2phase == 1) // 2 phase
{
   Info <<"in 2phase "<<endl;
	Ug.correctBoundaryConditions();
	forAll(mesh.boundary(),patchi)
	{
		if (isA< fixedValueFvPatchField<vector> >(Ug.boundaryField()[patchi]))
			{ phig.boundaryFieldRef()[patchi] = Ug.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi]; }
	} 
   #include "pSEqn.H"
}
else 
{
   Info <<"in simple flow "<<endl;
   #include "hEqn.H"
}
Uw = fvc::reconstruct(phiw); // This is a true darcy velocity (flux divided by surface even in unconfined)
