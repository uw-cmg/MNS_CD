#ifndef Constant_h
#define Constant_h

#define ZERO RCONST(0.0)        //zero
#define ONE RCONST(1.0)         //one
#define TWO RCONST(2.0)         //two

#define kb RCONST(8.617E-5)     //Boltzmann
#define pi RCONST(3.141592)     //pi

#define numPhase 2            	//Number of precipitating phases
#define numComp 3             	//Number of precipitating components
#define numClass 50000        	//Number of cluster classes/maximum cluster size considered
#define runs 70			//Number of loops to run
#define CutoffSize 65           //Cutoff size used for output
#define RadiusCalc radM2       	//Method to calculate mean radius of precipitate, radM1 or radM2
#define T0 ZERO                 //Initial time

#define numCalcPhase (numPhase*2)                       //number of calculating phases,including both home and heter nucleated phases
#define neq (numCalcPhase*numClass+numComp)             //number of ODEs

#define RTOL RCONST(1.0E-6)     //rel tolerance
#define ATOL RCONST(ZERO)       //abs tolerance

#endif
