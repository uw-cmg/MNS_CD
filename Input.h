#include <sundials/sundials_types.h> /* definition of type realtype */
#include "Constants.h"

#ifndef Input_h
#define Input_h
typedef struct InputMaterialType  {realtype aLat, aVol, C0[numComp],D[numComp],DFe,X[numPhase][numComp],aP[numPhase],cVol[numPhase],sig[numPhase],solPBar[numPhase];} *InputMaterial;
typedef struct InputConditionType {realtype Temp, Flux;} *InputCondition;
typedef struct InputPropertyType  {int HGSize, HGPhase; realtype Alpha, RsolP, ccs, Rflux, p_factor, DDP, DCB, SigmaDpa, DV, rv;} *InputProperty; //Cascade efficiency, v-i recombination radius,binding energy of di-vacancyies, dislocation densities 
void LoadInput( InputCondition ICond, InputMaterial IMaterial, InputProperty IProp);
#endif
