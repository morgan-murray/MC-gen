#ifndef PHYSICS_CONST_H
#define PHYSICS_CONST_H

#include <TMath.h>

const double ALPHA_EM = 0.00729735;
const double HC2 = 0.389268;
const double Mp = 0.938272;
const double Mele = 0.000510998;
const double Mmu = 0.105658;
const double Mpi = 0.139570;
const double PI = TMath::Pi();

const double Mp2 = pow(Mp,2);
const double Mele2 = pow(Mele,2);
const double Mmu2 = pow(Mmu,2);

const double coeff = pow(10,9)*HC2*pow(ALPHA_EM,3);

#endif
