#include <stdio.h>
#include <math.h>

#include "physics_const.h"
#include "physics_types.h"

inline double sq(double x) { return x*x;}

double calcBH1(struct _beam beam, struct _kinematics ks, 
	       struct _form_factors ff, double phig, double tau, double *dsBH);

void calcFFs(double del2, struct _form_factors * form_factors);

void calcCFFs(int gpdmod,int target, struct _kinematics kinematics, 
	      struct _form_factors form_factors, 
	      struct _compton_form_factors * cffs);

int bhdvcs(int xsecopt, int xsectyp, 
	   int gpdmod, int target,
	   struct _beam beam , struct _kinematics * kinematics, 
	   struct _exclusive_vectors * vectors,
	   double scattered_e_phi, double scattered_g_phi,
	   double * dsBH, double *dsDVCS, double *dsINT);

void calcBH(struct _beam beam, struct _kinematics kinematics, 
	    struct _form_factors form_factors, 
	    struct _exclusive_vectors * vectors,
	    double E_p, double E_Scatter, 
	    double E_gamma, double nu,
	    double qmod, double cos_t_gg, 
	    double tau, double lam1, double lam2, double produced_g_phi,
	    double *dsBH);

void calcDVCS(int gpdmod, int target, 
	      struct _beam beam, struct _kinematics * kinematics,
	      struct _form_factors form_factors,
	      double del2min, double del2max,
	      double lam1, double lam2, double tau, double phip,
	      double *dsDVCS, double *dsINT);

int readTab(int gpdmod,
	    double * re_h_u, double * re_h_d, 
	    double * re_del_u, double * re_del_d, 
	    double * re_e_u, double * re_e_d, 
	    double * re_pion_pole_u, double * re_pion_pole_d,
	    double * im_h_u, double * im_h_d, 
	    double * im_del_u, double * im_del_d,
	    double * im_e_u, double * im_e_d, 
	    double * im_pion_pole_u, double * im_pion_pole_d);
