#ifndef PHYSICS_TYPES_H
#define PHYSICS_TYPES_H

#include <stdio.h>
#include <TVector3.h>

struct _beam{

  int charge;
  double helicity;
  int particle; // 1 = electron, 2 = muon
  double energy;

};

struct _kinematics{

  double Q2;
  double xB;
  double del2;
  double y;
  double t;
};

struct _exclusive_vectors{

  TVector3 incidentLepton;
  TVector3 incidentTarget;
  TVector3 virtualPhoton;
  TVector3 outgoingLepton;
  TVector3 outgoingTarget;
  TVector3 outgoingPhoton;

};

struct _form_factors{
  
  double GE_p;
  double GE_n;
  double GM_p;
  double GM_n;

  double F1p;
  double F1n;
  double F2p;
  double F2n;
};

struct _compton_form_factors{

  double real_H;
  double imaginary_H;
  double real_Htilde;
  double imaginary_Htilde;
  
  double real_E;
  double imaginary_E;
  double real_Etilde;
  double imaginary_Etilde;
};

#endif
