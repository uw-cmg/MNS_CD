#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include "Input.h"
#include "Constants.h"
#include <sundials/sundials_types.h> /* definition of type realtype */


void LoadInput(InputCondition ICond, InputMaterial IMaterial, InputProperty IProp){

//*******************************************
//Input irradiation Condition
//*******************************************
	ICond->Temp=573;				//Temperature (K)
	ICond->Flux = 30E16; 				//Flux (m^-2s-1)


//*******************************************
//Input material property
//*******************************************
	IMaterial->aLat=2.87E-10;                       //Lattice constant (m)
	int Na=2;                                       //Number of atoms in cell, 2 for bcc, 4 for fcc
	IMaterial->aVol=pow(IMaterial->aLat,3)/Na;      //Atomic volume (m3)
	IMaterial->C0[0] = 0.0116; 			//Mn composition in alloy
        IMaterial->C0[1] = 0.0163; 			//Ni composition in alloy
        IMaterial->C0[2] = 0.0035; 			//Si composition in alloy
        IMaterial->D[0]= (1.49E-4)*exp(-234000.0/(8.314*ICond->Temp));  //Mn thermal diffusion coefficient, m2/s
        IMaterial->D[1] = (1.4E-4)*exp(-245765.0/(8.314*ICond->Temp));  //Ni thermal diffusion coefficient, m2/s
        IMaterial->D[2] = (0.78E-4)*exp(-231521.6/(8.314*ICond->Temp)); //Si thermal diffusion coefficient, m2/s
        IMaterial->DFe=(2.75E-3)*exp(-(253968.8)/(8.314*ICond->Temp));  //Fe thermal diffusion coefficient, m2/s
        IMaterial->X[0][0] = 6./29.; 			//Mn composition in T3/G (Mn6Ni16Si7) phase
        IMaterial->X[0][1] = 16./29.;			//Ni composition in T3/G (Mn6Ni16Si7) phase
        IMaterial->X[0][2] = 7./29.;			//Si composition in T3/G (Mn6Ni16Si7) phase
        IMaterial->X[1][0] = 1./3.;			//Mn composition in T6/Gamma2 (Mn(Ni,Si)2) phase
        IMaterial->X[1][1] = 0.5215;			//Ni composition in T6/Gamma2 (Mn(Ni,Si)2) phase
        IMaterial->X[1][2] = 0.1452;			//Si composition in T6/Gamma2 (Mn(Ni,Si)2) phase
        IMaterial->sig[0] = 0.185*6.241E18;		//T3/G phase interfacial energy is 0.185J/m2
        IMaterial->sig[1] = 0.175*6.241E18;             //T6/Gamma2 phase interfacial energy is 0.175J/m2
        IMaterial->cVol[0] = pow(1.093E-9,3)/116.;      //Atomic volume of T3/G phase, m^-3
        IMaterial->cVol[1] = pow(6.67E-10,3)/24.;       //Atomic volume of T6/Gamma2 phase,m^-3
        IMaterial->solPBar[0] = 2.45E-3;		//Equlibrium solute product of T3 phase at 573K
        IMaterial->solPBar[1] = 2.82E-3;		//Equilibrium solute product of T6 phase at 573K
//*******************************************
//Input other properties
//*******************************************

	IProp->HGSize=80;				//Heterogeneous nucleation size 
	IProp->HGPhase=1;				//Heterogeneous nucleation phase, 0 is T3 phase, 1 is T6 phase
	IProp->Alpha=4.8E-3;				//Cascade cluster production efficiency factor
	IProp->RsolP=2.4E-3;				//Reference solute product used in Eq. (5) in Sec. 2.2
	IProp->ccs=2E-28;				//Cascade cross section,m2
	IProp->Rflux=3E15;				//Reference flux, m^-2s^-1
        IProp->p_factor=0.2;				//p-factor used in Eq.(9) and (10) in Sec. 2.3
	IProp->DDP=2E14;				//Dislocation number density, m^-2
	IProp->DCB=0.4;             			//Cascade efficiency
	IProp->SigmaDpa=1.5E-25;   			//dpa cross section, m2
	IProp->DV=1E-4*exp(-1.3/(kb*ICond->Temp));	//Vacancy diffusion coefficient, m2/s
	IProp->rv=5.7e-10;				//SIA-vacancy recombination radius, m
}
