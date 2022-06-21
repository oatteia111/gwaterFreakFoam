#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"
int nthreads = 3;

//////////////////// find local dir
#include <unistd.h>
#define GetCurrentDir getcwd

std::string get_current_dir() {
   char buff[FILENAME_MAX]; //create string buffer to hold path
   GetCurrentDir( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
}
string cur_dir = get_current_dir();

///////////////// get some important values
/* ifstream f0(cur_dir/"phqfoam.txt");
std::vector<int> ph_data;
if (f0) {        
    int value;
    // read the elements in the file into a vector  
    while ( f0 >> value ) {
        ph_data.push_back(value);
    }
f0.close(); */

int nxyz,ph_ncomp,ph_gcomp,ph_nsolu,sunits;
////////////////////// start phreeqc
class my_phq
{
public:
	PhreeqcRM *PhreeqcRM_ptr;
#ifdef USE_MPI
	MPI_Comm rm_commxx;
#endif
	std::string DB;
	void setDB(std::string db){this->DB = db;}
	std::vector<int> data;
	void setData(std::vector<int> data){this->data = data;}
	std::string ChemFile;
	void setChemFile(std::string chfile){this->ChemFile = chfile;}
	std::vector<double> c;
	void setC(std::vector<double> c){this->c = c;}
	std::vector<double> g;
	void setG(std::vector<double> g){this->g = g;}
	std::vector<double> gm;
	void setGm(std::vector<double> gm){this->gm = gm;}
	std::vector<double> p;
	void setP(std::vector<double> p){this->p = p;}
	std::vector<double> spc;
	void setSpc(std::vector<double> spc){this->spc = spc;}
	std::vector<double> solu_conc;
	void setSconc(std::vector<double> solu_conc){this->solu_conc = solu_conc;}
	double tstep;
	void setTstep(double tstep){this->tstep = tstep;}
	std::vector<double> poro;
	void setPoro(std::vector<double> poro){this->poro = poro;}
	std::vector<double> wsat;
	void setWsat(std::vector<double> wsat){this->wsat = wsat;}
	std::vector<double> gvol;
	void setGvol(std::vector<double> gvol){this->gvol = gvol;}
	int iGwater=-1;
	void setNwat(int){this->iGwater = iGwater;}

	//std::vector<std::string> complist;
    //void setCompList(std::vector<std::string>) {this->complist = complist;}
	//this->complist = complist;
	void init()
	{
		nxyz = this->data[0];ph_ncomp = this->data[1];ph_gcomp = this->data[2];sunits = this->data[4];std::cerr << " ph nxyz "<< nxyz <<"\n";
		int i,j;
		std::vector<double> c(nxyz*(ph_ncomp+ph_gcomp),0.);
		this->c = c;

		PhreeqcRM *phreeqc_rm = new PhreeqcRM (nxyz, nthreads);
		this->PhreeqcRM_ptr = phreeqc_rm;
		IRM_RESULT status;
		status = phreeqc_rm->SetErrorHandlerMode(1);
		status = phreeqc_rm->SetComponentH2O(true); // to get water,excess H an dexcess O transported
		status = phreeqc_rm->SetRebalanceFraction(0.5);
		status = phreeqc_rm->SetRebalanceByCell(true);
		phreeqc_rm->UseSolutionDensityVolume(false);
		phreeqc_rm->SetPartitionUZSolids(false);
		status = phreeqc_rm->SetUnitsSolution(2);           // 1, mg/L; 2, mol/L; 3, kg/kgs
		status = phreeqc_rm->SetUnitsPPassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsExchange(sunits);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSurface(sunits);            // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsGasPhase(sunits);           // 0, mol/L RV; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsSSassemblage(sunits);       // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
		status = phreeqc_rm->SetUnitsKinetics(sunits);           // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
				
		std::vector<double> rv;
		rv.resize(nxyz, 1.);
		status = phreeqc_rm->SetRepresentativeVolume(rv);   // RV is one dm3 of medium
		Info<<"poro size "<<this->poro.size()<<endl;
		/*std::vector<double> por;
		por.resize(nxyz, 0.25);
		status = phreeqc_rm->SetPorosity(por);*/
		status = phreeqc_rm->SetPorosity(this->poro);
		Info << "end setUnits... poro 10 350 " << this->poro[10] << " " << this->poro[350] << endl;
		
		// Open for writing file
		status = phreeqc_rm->SetFilePrefix("phq_out");
		phreeqc_rm->OpenFiles();
		//Info << "Outfile opened " << endl;
		//set mask
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 1);
		//for (int i=0;i<nxyz;i++) {print_chemistry_mask[i]=1;}
		status = phreeqc_rm->SetPrintChemistryMask(print_chemistry_mask);

		//load database and chemistry files and make the first run of the chem file
		status = phreeqc_rm->SetPrintChemistryOn(false, true, false); // workers, initial_phreeqc, utility
		status = phreeqc_rm->LoadDatabase(this->DB);
	Info << "dbase opened " << endl;
		bool workers = true;             // Worker instances do the reaction calculations for transport
		bool initial_phreeqc = true;     // InitialPhreeqc instance accumulates initial and boundary conditions
		bool utility = true;             // Utility instance is available for processing
		status = phreeqc_rm->RunFile(workers, initial_phreeqc,utility, this->ChemFile);

		// Clear contents of workers and utility
		initial_phreeqc = false;
		std::string input = "DELETE; -all";
		status = phreeqc_rm->RunString(workers, initial_phreeqc, utility, input.c_str());
		
		int ph_ncomp = phreeqc_rm->FindComponents(); // NECESSARY don't remove

		// Set array of initial conditions by reading the solution... numbers in data
		std::vector<int> ic1, ic2;
		ic1.resize(nxyz*7, -1);
		//std::vector<double> f1;
		int icnt = 4; // to consider the first 5 numbers : nxyz,ncomp,gcomp,nsolu,units (++ in the first line of the for makes 5 numbers)
		for (i = 0; i < 7; i++) // Solution,Equilibrium phases,Exchange,Surface,Gas phase,Solid solutions,Kinetics
		{
			icnt++;
			int test = this->data[icnt];//Info<<"index "<<i<<" " <<test<<endl;
			if (test == 0)  // if 0 just one value to read
			{icnt++; for (int j=0;j<nxyz;j++) {ic1[i*nxyz+j] = this->data[icnt];
				if (j==369) {Info << i << " cell "<< j << " " << this->data[icnt]<< endl;}
			}}
			else 		
			{for (int j=0;j<nxyz;j++) {icnt++;ic1[i*nxyz+j] = this->data[icnt];
				if (j==369) {Info << i << " cell "<< j << " " << this->data[icnt]<< endl;}
			}}
		}
		status = phreeqc_rm->InitialPhreeqc2Module(ic1);
		
		// Set solutions composition at fixed cells HERE only TWO solutions (shall be modified)
		//std::vector<double> solu_conc;
		//std::vector<int> cf1 ={0,-1};
		//for (i= 0;i<2;i++) {cf1[this->cfix[2+i*4+2]] = this->cfix[2+i*4+3];Info<<" cfix "<<i<<" "<< this->cfix[2+i*4+2] << endl;}

		// run initial conditions
		int ncell = phreeqc_rm->GetChemistryCellCount();
		std::cerr << "ncell "<< ncell << " ph_ncomp "<< ph_ncomp<<" starting init phreeqc \n";
		//status = phreeqc_rm->SetGasPhaseVolume(this->gvol);
		//status = phreeqc_rm->SetPressure(this->p);
		//phreeqc_rm->SetPrintChemistryOn(true,false,false);
		status = phreeqc_rm->RunCells();
		std::cerr << " cells run \n";
		status = phreeqc_rm->GetConcentrations(this->c);
		status = phreeqc_rm->GetGasCompPressures(this->g);
		for (i=0;i<ph_gcomp;i++){for (j=0;j<5;j++) {Info <<i<<" "<<j << " g "<< this->g[i*nxyz+j]<< endl;} }
		//status = phreeqc_rm->InitialPhreeqc2Concentrations(this->solu_conc, cf1);		
		const std::vector<std::string> &gcomp = phreeqc_rm->GetGasComponents(); //Info<<"phq gcomp0 "<<gcomp[0]<<endl;
		for (i=0;i<gcomp.size();i++) {
			if (gcomp[i]=="H2O(g)") {this->iGwater = i;}
			}


		/* Set boundary condition
		std::vector<double> bc_conc;
		std::vector<double> bc_f1; //bc_conc, 
		std::vector<int> bc1, bc2;
		int nbound = 2;
		bc1.resize(nbound, 0);                      // solution 0 from Initial IPhreeqc instance to boundary condition
		bc2.resize(nbound, 1);                     // no bc2 solution for mixing
		bc_f1.resize(nbound, 1.0);                  // mixing fraction for bc1
		status = phreeqc_rm->InitialPhreeqc2Concentrations(this->solu_conc, bc1, bc2, bc_f1); */
		
		
		//printing to file
		status = phreeqc_rm->CloseFiles();
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, true); // workers, initial_phreeqc, utility

	}
	void run()
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
		/*
		status = phreeqc_rm->SetErrorHandlerMode(1);
		status = phreeqc_rm->SetFilePrefix("phq_run");
		phreeqc_rm->OpenFiles();
		std::vector<int> print_chemistry_mask;
		print_chemistry_mask.resize(nxyz, 0);print_chemistry_mask[0]=1;print_chemistry_mask[50]=1;
		status = phreeqc_rm->SetPrintChemistryMask(print_chemistry_mask); */
		//status = phreeqc_rm->SetPrintChemistryOn(true, false, false); // workers, initial_phreeqc, utility
		const std::vector<std::string> &components = phreeqc_rm->GetComponents();
		const std::vector<std::string> &gcomp = phreeqc_rm->GetGasComponents(); //Info<<"phq gcomp0 "<<gcomp[0]<<endl;
		status = phreeqc_rm->SetSaturation(this->wsat);
		status = phreeqc_rm->SetConcentrations(this->c);
		if (gcomp.size()>0) {
			status = phreeqc_rm->SetGasPhaseVolume(this->gvol);
			status = phreeqc_rm->SetGasCompMoles(this->gm);
			}
	Info << "in phq poro "<<this->poro[2]<<" sw 2 "<< this->wsat[2] << " gvol "<< this->gvol[2] <<" g "<<this->g[2]<< endl;
	for (int i=0; i<ph_ncomp;i++){Info << "c 0 "<<i<<" "<< this->c[i*nxyz] << endl;}
	for (int i=0; i<ph_gcomp;i++){Info << "gm 0 "<<i<<" "<< this->gm[i*nxyz] << endl;}
		if (gcomp.size()>0) {
			status = phreeqc_rm->GetSaturation(this->wsat);
			status = phreeqc_rm->GetGasCompPressures(this->g);
			status = phreeqc_rm->GetGasCompMoles(this->gm);
			status = phreeqc_rm->GetGasPhaseVolume(this->gvol);
			}
	for (int i=0; i<ph_gcomp;i++){Info << "g 0 "<<i<<" "<< this->g[i*nxyz] << endl;}
	for (int i=0; i<ph_gcomp;i++){Info << "gm 0 "<<i<<" "<< this->gm[i*nxyz] << endl;}
	Info << "wsat 0 "<< this->wsat[0] << endl;
	Info << "gvol 0 "<< this->gvol[0] << endl;
	//status = phreeqc_rm->SetPressure(this->p);
		status = phreeqc_rm->SetTimeStep(this->tstep);
		status = phreeqc_rm->RunCells();
		status = phreeqc_rm->GetSaturation(this->wsat);
		status = phreeqc_rm->GetConcentrations(this->c);
		if (gcomp.size()>0) {
			status = phreeqc_rm->GetGasCompPressures(this->g);
			status = phreeqc_rm->GetGasCompMoles(this->gm);
			status = phreeqc_rm->GetGasPhaseVolume(this->gvol);
			}
		//writes to the dump file
		bool dump_on = true;
		bool append = false;
		status = phreeqc_rm->SetDumpFileName("phreqc.dmp");
		status = phreeqc_rm->DumpModule(dump_on, append);
	Info << "in phq poro "<<this->poro[2]<<" sw 2 "<< this->wsat[2] << " gvol "<< this->gvol[2] <<" g "<<this->g[2]<< endl;
	//for (int i=0; i<ph_ncomp;i++){Info << "c middle "<<i<<" "<< this->c[i*nxyz+nxyz/2] << endl;}
	}
	void getSelOutput()
	{
		IRM_RESULT status;
		PhreeqcRM* phreeqc_rm = this->PhreeqcRM_ptr;
		int n_user = phreeqc_rm->GetNthSelectedOutputUserNumber(0);
		status = phreeqc_rm->SetCurrentSelectedOutputUserNumber(n_user);
		int ncol = phreeqc_rm->GetSelectedOutputColumnCount();
		/* std::string heading;
		for (int i = 0; i < ncol; i++)
			{
				status = phreeqc_rm->GetSelectedOutputHeading(i,heading);
				if (heading=="g_H2O(g)") {this->nwat=i-6;Info<<heading<<" nb wat "<<this->nwat<<endl;}
			} */
		nxyz = this->data[0];
		this->spc.resize((ncol-6)*nxyz,0);
		std::vector<double> so;
		status = phreeqc_rm->GetSelectedOutput(so);
		//std::vector<std::string> headings;
		//headings.resize(ncol);
		for (int i = 0; i < nxyz; i++)
		{
			for (int j = 6; j < ncol; j++)
			{
				this->spc[(j-6)*nxyz + i] = so[j*nxyz + i];
				if (i==10) {Info<< "selout sp "<< j << " : " << so[j*nxyz + i] << endl;}
			}
		}
	}
};
