//Code developed by Huibin Ke from UW-Madison for the evolution of Mn-Ni-Si precipitates in Reactor Pressure Vessel stees.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <cvode/cvode.h>             /* main integrator header file */
#include <cvode/cvode_band.h>       /* band solver header */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_types.h> /* definition of realtype */
#include <sundials/sundials_math.h>  /* contains the macros ABS, SQR, and EXP */

#include <iostream>

#include "Input.h"          /*Input parameter header file*/
#include "Constants.h"      /*Constants header file*/


InputCondition ICond;      /*Defined in Input.h, including irradiation conditions*/
InputMaterial IMaterial;   /*Defined in Input.h, including material information*/
InputProperty IProp;       /*Defined in Input.h, including all other parameters used in model*/

using namespace std;

/*function defs*/

typedef struct {
  realtype size[numClass],radClust[numPhase][numClass],beta[numPhase][numClass],delG[numPhase][numClass];
} *UserData;

static void loadData(UserData& data);
static void getInitVals(realtype y0[neq]);
static void initParams();
static void GetRED(realtype D[numComp], realtype Flux);
static void printYVector(N_Vector y);
static void getOutput(N_Vector y, realtype radM1[numCalcPhase], realtype radM2[numCalcPhase],realtype rhoC[numCalcPhase]);
static void getbeta(realtype size[numClass], realtype beta[numPhase][numClass]);
static void getSize(realtype size[numClass]);
static void getRadClust(realtype size[numClass], realtype radClust[numPhase][numClass]);
static void getDelG(realtype size[numClass], realtype radClust[numPhase][numClass], realtype delG[numPhase][numClass]);
static void getFlux(UserData data, N_Vector y, realtype J[numCalcPhase][numClass+1]);
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void int_to_string(int i, string &a, int base);

/*global problem parameters*/
realtype D[numComp],aP[numPhase],nu[numPhase][numComp]; /*Diffusion coefficient,effective precipitate lattice constant, square of precipitate composition*/
realtype rhoC[numCalcPhase],radM1[numCalcPhase],radM2[numCalcPhase];/*precipitate number density, two kinds of mean radius (see readme for more detail)*/
realtype Flux,solProd[numPhase]; /*Irradiation flux, solute product */

int main(){

	realtype y0data[neq],t;
	realtype tout = 1E0; //This is the solution time of each run.  The solver will solve for y0(t) on the range (0,tout)
	double ts=0.0; //This is the start time for each run, tf+tout will be the final time ater each run

	/*The following block setup the information for the solver and initial all the parameters*/
	N_Vector y0 = NULL;
        ICond=NULL;
	IMaterial=NULL;
	IProp=NULL;
	UserData data = NULL;
	void *cvode_mem = NULL;
	int flag,itask,mxsteps;
	mxsteps = 2000000;
        long int mu = 3;
	long int ml = 3;
	realtype *yd;
	ICond=(InputCondition)malloc(sizeof *ICond);         //Allocate irradiation condition memory
	IMaterial=(InputMaterial)malloc(sizeof *IMaterial);  //Allocate material properties memory
	IProp=(InputProperty)malloc(sizeof *IProp);          //Allocate other properties memory
	LoadInput(ICond,IMaterial,IProp);	             //Load all the input parameters
	initParams();
        getInitVals(y0data);
	y0 = N_VMake_Serial(neq,y0data);


	/*Write the output file*/
	ofstream O_file;
        O_file.open("Output");
        O_file << "Run\tCalcTime(s)\tTime(s)\tFluence(n/m2s)\t";
	for(int p=0; p<numPhase; p++){
		string phaseStr;
		int_to_string(p+1,phaseStr,10);
		cout<<p<<" "<<p+1<<" "<<phaseStr<<endl;
        	O_file<<"Mean_Radius_of_Phase_"+phaseStr+"_Homo(m)\tNumber_Density_of_Phase_"+phaseStr+"_Homo(1/m3)\t";
	}
	for(int p=0; p<numPhase; p++){
		string phaseStr;
		int_to_string(p+1,phaseStr,10);
                O_file<<"Mean_Radius_of_Phase_"+phaseStr+"_Heter(m)\tNumber_Density_of_Phase_"+phaseStr+"_Heter(1/m3)\t";
        }
        O_file <<"Mn\tNi\tSi"<<endl;


	/*Solving all the equations, each run gives output at time ts+tout*/
	for(int i=0; i<runs; i++){
     		cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);
		data = (UserData) malloc(sizeof *data);
		loadData(data);
		flag = CVodeSetUserData(cvode_mem, data);
		flag = CVodeInit(cvode_mem, f, T0, y0);
		flag = CVodeSStolerances(cvode_mem,RTOL,ATOL);
		flag = CVBand(cvode_mem, neq, mu, ml );
		flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps);
		time_t tik,tok;  /*tik is the time before solving equations, tok is the time after solving equations, tok-tik will be the time used in each run*/
		time(&tik);
		flag = CVode(cvode_mem, tout, y0, &t, CV_NORMAL);
		time(&tok);
		getOutput(y0, radM1, radM2, rhoC); /*Obtain the output after each run*/
		yd = NV_DATA_S(y0);
		/*The following few lines write the output into file Output*/
		O_file << i << "	" <<difftime(tok,tik) << "	" << ts+tout << "	"<< (ts+tout)*Flux<<"	";
		for(int p=0; p<numCalcPhase; p++){
			O_file << RadiusCalc[p] << "	" << rhoC[p] << "	";
		}
		O_file << yd[neq-3] << "	" << yd[neq-2] << "	" << yd[neq-1];
		O_file << endl;
		printYVector(y0);  /*Write the particle size distribution of each phase into a seperate file*/
		CVodeFree(&cvode_mem);
		ts = ts + tout;     /*Calculate the start time of next run*/
		tout = pow(10,i/9); /*Define the solution time for next run*/
	}
	N_VDestroy_Serial(y0);
return 0;
}

/*************************************************

This function initialize parameters used in model

*************************************************/
static void initParams(){
	Flux = ICond->Flux; /*Irradiatio flux*/
	for(int i=0; i<numComp;i++){
		D[i]=IMaterial->D[i];    /*Thermal diffusion coefficients*/
	}
	GetRED(D,Flux);	 /*Calculate radiation enhanced diffusion coefficients*/
        for(int p=0; p<numPhase; p++){
		aP[p] = pow((3*IMaterial->cVol[p])/(4*pi),1./3.);   /*Effective atomic radius in precipitate*/
	}
	for(int p=0; p<numPhase; p++){
		for(int c=0; c<numComp; c++){
			nu[p][c] = pow(IMaterial->X[p][c],2);  /*Square of precipitate composition*/
		}
	}

return;
}

/*****************************************************************

This function gives the radiation enhanced diffusion coefficients
Described in SI Sec. C

*****************************************************************/
static void GetRED(realtype D[numComp], realtype Flux){
        realtype Eta,Gs,Cvr;

	/*When flux is higher than reference flux, use p-scaling (Eq. SI-18) to calculate gs*/
	if (Flux>IProp->Rflux){
		Eta = 16*pi*IProp->rv*IProp->DCB*IProp->SigmaDpa*IProp->Rflux/IMaterial->aVol/IProp->DV/pow(IProp->DDP,2);
		Gs = 2.0/Eta*(pow(1+Eta,0.5)-1.0)*pow((IProp->Rflux/Flux),IProp->p_factor);
	} /*When flux is lower than reference flux, use Eq. SI-19 and SI-20 to calculate gs*/
	else {
		Eta = 16*pi*IProp->rv*IProp->DCB*IProp->SigmaDpa*Flux/IMaterial->aVol/IProp->DV/pow(IProp->DDP,2);
		Gs = 2.0/Eta*(pow(1+Eta,0.5)-1.0);
	}
	Cvr = IProp->DCB*Flux*IProp->SigmaDpa*Gs/(IProp->DV*IProp->DDP); /*Calculate vacancy concentration under irradiation, Eq. SI-17*/
        for(int i=0;i<numComp; i++){
		D[i]=D[i]+IProp->DV*Cvr*D[i]/IMaterial->DFe; /*Radiation enhanced diffusion coefficients, Eq. SI-16*/
	}
}


/*****************************************************************

This function loads data that is as a function of cluster size
All values are calcuate via next few functions.

*****************************************************************/
static void loadData(UserData &data){
	realtype size[numClass],radClust[numPhase][numClass],beta[numPhase][numClass],delG[numPhase][numClass];
	getSize(size);
        getRadClust(size,radClust);
        getbeta(size,beta);
        getDelG(size,radClust,delG);

	for(int i=0; i<numClass; i++){
		data->size[i] = size[i];
	}

	for(int p=0; p<numPhase; p++){
		for(int i=0; i<numClass; i++){
			data->radClust[p][i] = radClust[p][i];
			data->beta[p][i] = beta[p][i];
			data->delG[p][i] = delG[p][i];
		}
	}

return;
}

/*****************************************************************

This function calculates part of the absorption rate

*****************************************************************/
static void getbeta(realtype size[numClass], realtype beta[numPhase][numClass]){
	for(int p=0; p<numPhase; p++){
		for(int i=0; i<numClass; i++){
			beta[p][i] = (4*pi*aP[p]*pow(size[i],(1./3.)))/IMaterial->aVol;
		}
	}

return;
}

/*****************************************************************

This function calculates number of atoms (size) in each cluster

*****************************************************************/

static void getSize(realtype size[numClass]){
	for(int i=0; i<numClass; i++){
		size[i] = double(i+1);
	}
return;
}

/*****************************************************************

This function calculates radius of each cluster

*****************************************************************/

static void getRadClust(realtype size[numClass], realtype radClust[numPhase][numClass]){

	for(int p=0; p<numPhase; p++){
		for(int i=0; i<numClass; i++){
			radClust[p][i] = pow(3*IMaterial->cVol[p]*size[i]/(4*pi),1./3.);
		}
	}
return;
}


/*********************************************************************************

This function calculates difference of interfacial energy between adjacent clusters

**********************************************************************************/

static void getDelG(realtype size[numClass], realtype radClust[numPhase][numClass], realtype delG[numPhase][numClass]){
	realtype delta;

	for(int p=0; p<numPhase; p++){
		for(int i=1; i<numClass; i++){
			delta = -(
				(IMaterial->sig[p]*4*pi*pow(radClust[p][i-1],TWO))-
				(IMaterial->sig[p]*4*pi*pow(radClust[p][i],TWO))
				 );
			delG[p][i] = exp(delta/(kb*ICond->Temp));
		}
	}
return;
}

/*****************************************************************

This function gives the inital values of number density of each cluster
1E-30 for clusters >= 2 atoms
monomer concentration is calculated based on Eq. SI-15

*****************************************************************/

static void getInitVals(realtype y0[neq]){
	realtype base = 1E-30;
	for(int i=0; i<neq; i++){
		y0[i] = base;
	}
	for(int c=0; c<numComp; c++){
		y0[neq-numComp+c] = IMaterial->C0[c];  /*Concentration of solute in matrix*/
	}
	for(int p=0; p<numPhase; p++){
		solProd[p] = 1;
		for(int c=0; c<numComp; c++){
			solProd[p] = solProd[p]*pow(IMaterial->C0[c],IMaterial->X[p][c]);
		}
		y0[p*numClass] = solProd[p];   /*Effective monomer concentration, based on Eq. SI-15*/
	}
return;
}

/**************************************************************************

This function calculates flux (Jn->n+1) between adjacent clusters, Eq. SI-2

****************************************************************************/

static void getFlux(UserData data, N_Vector y, realtype J[numCalcPhase][numClass+1]){
	realtype *yd;
	yd = NV_DATA_S(y);
	realtype solP[numCalcPhase],wp[numComp],sumwp,wpEff;
        int pref;   /*pref is used to refer to the real phase for both homo and heter nucleated phases, pref=0 is T3 phase, pref=1 is T6 phase*/
	for(int p=0; p<numPhase; p++){
		solP[p] = 1;
		for(int c=0; c<numComp; c++){
			solP[p] = solP[p]*pow(yd[neq-numComp+c],IMaterial->X[p][c]); /*Calculate solute product*/
		}
	}
  	for(int p=numPhase;p<numCalcPhase;p++){
		solP[p]=1E-30;                      /*solute product for heterogeneous nucleation phase, used to calcuate effective monomer concentration*/
	}
       	for(int p = 0; p<numCalcPhase; p++){
		pref = p%numPhase;
		J[p][0] = ZERO;
		yd[p*numClass] = solP[p];          /*Effective monomer concentration*/
		sumwp = 0;
                for(int c=0; c<numComp; c++){
                        wp[c] = yd[neq-numComp+c]*D[c]*data->beta[pref][0];
                        sumwp = sumwp+(nu[pref][c]/wp[c]);
                }
                wpEff = 1.0/sumwp;     /*Absorption rate for momoner to dimer*/
                J[p][1] = wpEff*(yd[p*numClass]/numPhase-(IMaterial->solPBar[pref]/solP[pref])*data->delG[pref][1]*yd[p*numClass+1]);	/*flux from monomer to dimer*/
		for(int i=2; i<numClass; i++){
			sumwp = 0;
			for(int c=0; c<numComp; c++){
				wp[c] = yd[neq-numComp+c]*D[c]*data->beta[pref][i-1];
				sumwp = sumwp+(nu[pref][c]/wp[c]);
			}
			wpEff = 1.0/sumwp;
			J[p][i] = wpEff*(yd[p*numClass+i-1]-(IMaterial->solPBar[pref]/solP[pref])*data->delG[pref][i]*yd[p*numClass+i]); /*flux from size i to size i+1*/
		}
		J[p][numClass] = ZERO;
	}
return;
}

/**********************************************************************************

This function calculates Eq. SI-1, the change of concentration of each cluster size

***********************************************************************************/

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data){
	realtype *yd,*ydotd;
	yd = NV_DATA_S(y);
	ydotd = NV_DATA_S(ydot);
	UserData data;

	realtype size[numClass],radClust[numPhase][numClass],beta[numPhase][numClass],delG[numPhase][numClass],J[numCalcPhase][numClass+1];
	int pref;
	realtype solP[numPhase];
	realtype GR;

	data = (UserData) user_data;

	getFlux(data,y,J);

	realtype sumNdot[numCalcPhase],Cdot;
	for(int p=0; p<numCalcPhase; p++){
		sumNdot[p] = ZERO;
		ydotd[p*numClass] = ZERO;
		for(int i=1; i<numClass; i++){
			ydotd[p*numClass+i] = J[p][i]-J[p][i+1];     /*Eq. 1 in Sec. 2.1 without Rhet term*/
			sumNdot[p] = sumNdot[p]+ydotd[p*numClass+i]*data->size[i]; /*Monomer consumed*/
		}
	}
        for(int p=numPhase;p<numCalcPhase;p++){
		pref = p%numPhase;
		solP[pref] = 1;
                for(int c=0; c<numComp; c++){
                        solP[pref] = solP[pref]*pow(yd[neq-numComp+c],IMaterial->X[pref][c]); /*solute product*/
                }
	}

	/*The next three lines is the generation of clusters in cascade damage*/
	GR = IProp->Alpha*Flux*IProp->ccs*solP[IProp->HGPhase]/IProp->RsolP;  /*Generation rate of clusters, Eq. 5 in Sec. 2.2*/
        ydotd[(IProp->HGPhase+numPhase)*numClass+IProp->HGSize-1]+=GR;        /*Add Rhet term to Eq. (1) in Sec. 2.1*/
        sumNdot[(IProp->HGPhase+numPhase)] += GR*data->size[IProp->HGSize-1]; /*Calculate the monomers consumed in heterogeneous nucleation*/

	for(int c=0; c<numComp; c++){
		Cdot = 0;
		for(int p=0; p<numCalcPhase; p++){
			pref = p%numPhase;
			Cdot = Cdot+IMaterial->X[pref][c]*sumNdot[p];
		}
		ydotd[neq-numComp+c] = -Cdot; /*change of solute/monomer in matrix*/
	}

return(0);
}


/*****************************************************************

This function calculates mean cluster radius and cluster density

*****************************************************************/

static void getOutput(N_Vector y, realtype radM1[numCalcPhase], realtype radM2[numCalcPhase],realtype rhoC[numCalcPhase]){
        realtype *yd,size[numClass],radClust[numPhase][numClass],numC,numCxSize1,numCxSize2;
        yd = NV_DATA_S(y);
	int pref;
        getSize(size);
        getRadClust(size,radClust);
	for(int p=0; p<numCalcPhase; p++){
		pref = p%numPhase;
        	numC = 0;            /*Concertration of clusters in unit per lattice site*/
        	numCxSize1 = 0;
                numCxSize2 = 0;
        	radM1[p] = 0;        /*Mean radius of clusters calculated via sum(r_i*N_i)/sum(N_i)*/  /*refer to readme for details regarding these two different ways of calculating mean radius*/
 		radM2[p] = 0;        /*Mean radius of clusters calculated via equilivant radius of a cluster with size sum(n_i*N_i)/sum(N_i)*/
        	rhoC[p] = 0;
        	for(int i=CutoffSize; i<numClass; i++){   /*Count clusters larger than CutoffSize for output*/
                	numC = numC+yd[p*numClass+i];
                	numCxSize1 = numCxSize1+yd[p*numClass+i]*radClust[pref][i];
                        numCxSize2 = numCxSize2+yd[p*numClass+i]*size[i];
        	}
                radM1[p] = numCxSize1/numC;
        	radM2[p] = pow((numCxSize2/numC*IMaterial->cVol[pref])/((4./3.)*pi),(1./3.));
        	rhoC[p] = numC/IMaterial->aVol; /*Calculate number density of clusters in unit of per m^3*/
	}
return;
}

/**********************************************************************************************

This function prints cluster size distribution in the file Profile for the final solution time.

***********************************************************************************************/
static void printYVector(N_Vector y){
        realtype *yd,size[numClass],radClust[numPhase][numClass];
        ofstream P_file;
        yd = NV_DATA_S(y);
	int pref;

	for(int p=0; p<numCalcPhase; p++){
		pref = p%numPhase;
		string profStr = "Profile_";
		string phaseStr;
		int_to_string(p,phaseStr,10);
		profStr.append(phaseStr);
        	P_file.open(profStr.c_str());
		P_file << "cluster size (# atoms)	cluster radius (m)	cluster density (1/m3)" << endl;
        	getSize(size);
        	getRadClust(size,radClust);
        	for(int i=0; i<numClass; i++){
                	P_file << size[i] << "	" << radClust[pref][i] << "	" << yd[p*numClass+i]/IMaterial->aVol << endl;
        	}
        	P_file.close();
	}
return;
}

void int_to_string(int i, string &a, int base){
  int ii=i;
  string aa;
  int remain=ii%base;

  if(ii==0){
    a.push_back(ii+48);
    return;
  }

  while(ii > 0){
    aa.push_back(ii%base+48);
    ii=(ii-remain)/base;
    remain=ii%base;
  }
  for(ii=aa.size()-1; ii>=0; ii--){
    a.push_back(aa[ii]);
  }
  return;
}
