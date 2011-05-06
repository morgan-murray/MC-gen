#include "bhdvcs.h"

void calcFFs(double del2, struct _form_factors * form_factors){
  
  const double Mv = 0.84;
  const double kp = 1.79; // Magnetic moment proton
  const double kn = -1.91; // Magnetic moment neutron
   
  const double dipol = 1.0/pow((1.0 - del2/pow(Mv,2)),2);
  
  form_factors->GE_p = dipol;
  form_factors->GE_n = 0;
  form_factors->GM_p = (1.0 +kp)*dipol;
  form_factors->GM_n = kn*dipol;
  
  const double delm = del2/(4*Mp2);
  const double denom = (1.0 - delm);

  form_factors->F1p = (form_factors->GE_p - delm*form_factors->GM_p)/denom;
  form_factors->F1n = (form_factors->GE_n - delm*form_factors->GM_n)/denom;
  form_factors->F2p = (form_factors->GM_p - form_factors->GE_p)/denom;
  form_factors->F2n = (form_factors->GM_n - form_factors->GE_n)/denom;

  return;
}
  
// Calculate the CFFs based on Korotkov's model in hep-ph/0108077v1
void calcCFFs(int gpdmod, int target, struct _kinematics kinematics, 
	      struct _form_factors form_factors, struct _compton_form_factors * cffs){

  // Values to be read from gpd.dat file
  double re_h_u[51], re_h_d[51], re_del_u[51],re_del_d[51]; 
  double re_e_u[51], re_e_d[51], re_pion_pole_u[51], re_pion_pole_d[51];
  double im_h_u[51], im_h_d[51], im_del_u[51], im_del_d[51];
  double im_e_u[51], im_e_d[51], im_pion_pole_u[51], im_pion_pole_d[51];


  readTab(gpdmod,
	  re_h_u, re_h_d ,re_del_u, re_del_d,
	  re_e_u, re_e_d, re_pion_pole_u, re_pion_pole_d,
	  im_h_u, im_h_d ,im_del_u, im_del_d,
	  im_e_u, im_e_d, im_pion_pole_u, im_pion_pole_d);

  double F1u = 2.0*form_factors.F1p + form_factors.F1n;
  double F1d = form_factors.F1p + 2.0*form_factors.F1n;
  double F2u = 2.0*form_factors.F2p + form_factors.F2n;
  double F2d = form_factors.F2p + 2.0*form_factors.F2n;

  double skewedness = kinematics.xB/(2 - kinematics.xB);

  // pdf values are in a log scale
  double min_skew = 0.01;
  double max_skew = 1.0;
  double log_skewedness = log10(skewedness);
  double log_min_skew = log10(min_skew);
  double log_max_skew = log10(max_skew);
  
  double skewedness_bins = (log_max_skew - log_min_skew)/51.0;
  
  // Interpolate between bins
  double select_bin = (log_skewedness - log_min_skew)/skewedness_bins;
  int low_bin = (int) floor(select_bin);
  double bin_interpol1 = select_bin - low_bin;
  double bin_interpol2 = 1 - bin_interpol1;

  // Get the pdf values at the interpolated value between bin edges
  double interpol_h_u_re = re_h_u[low_bin]*bin_interpol2 + re_h_u[low_bin+1]*bin_interpol1;
  double interpol_h_d_re = re_h_d[low_bin]*bin_interpol2 + re_h_d[low_bin+1]*bin_interpol1;

  double interpol_del_u_re = re_del_u[low_bin]*bin_interpol2 + re_del_u[low_bin+1]*bin_interpol1;
  double interpol_del_d_re = re_del_d[low_bin]*bin_interpol2 + re_del_d[low_bin+1]*bin_interpol1;

  double interpol_e_u_re = re_e_u[low_bin]*bin_interpol2 + re_e_u[low_bin+1]*bin_interpol1;
  double interpol_e_d_re = re_e_d[low_bin]*bin_interpol2 + re_e_d[low_bin+1]*bin_interpol1;

  double interpol_pion_pole_u_re = re_pion_pole_u[low_bin]*bin_interpol2 + re_pion_pole_u[low_bin+1]*bin_interpol1;
  double interpol_pion_pole_d_re = re_pion_pole_d[low_bin]*bin_interpol2 + re_pion_pole_d[low_bin+1]*bin_interpol1;

  double interpol_h_u_im = im_h_u[low_bin]*bin_interpol2 + im_h_u[low_bin+1]*bin_interpol1;
  double interpol_h_d_im = im_h_d[low_bin]*bin_interpol2 + im_h_d[low_bin+1]*bin_interpol1;

  double interpol_del_u_im = im_del_u[low_bin]*bin_interpol2 + im_del_u[low_bin+1]*bin_interpol1;
  double interpol_del_d_im = im_del_d[low_bin]*bin_interpol2 + im_del_d[low_bin+1]*bin_interpol1;

  double interpol_e_u_im = im_e_u[low_bin]*bin_interpol2 + im_e_u[low_bin+1]*bin_interpol1;
  double interpol_e_d_im = im_e_d[low_bin]*bin_interpol2 + im_e_d[low_bin+1]*bin_interpol1;

  // Now use axial-vector form factors to calculate the correct values for the cff
  double gA = 1.267/pow((1-kinematics.del2/0.84),2);
  double gA0 = 0.6*gA;
  double gAu = 0.5*(gA + gA0)/(0.8*1.267);
  double gAd = 0.5*(gA0 - gA)/(-0.2*1.267);

  // pion pole function
  double hA = 4*Mp2*1.267/(pow(Mpi,2)-kinematics.del2);

  // fill the cffs with the appropriate value
  if(target==1){

    cffs->real_H             = (F1u*interpol_h_u_re*2.0 + F1d*interpol_h_d_re)/9.0;
    cffs->imaginary_H        = (F1u*interpol_h_u_im*2.0 + F1d*interpol_h_d_im)/9.0;
    
    cffs->real_Htilde        = (4*gAu*interpol_del_u_re + gAd*interpol_del_d_re)/9.0;
    cffs->imaginary_Htilde   = (4*gAu*interpol_del_u_im + gAd*interpol_del_d_im)/9.0;
    
    cffs->real_E             = (F2u*interpol_e_u_re*2.0 + F2d*interpol_e_d_re)/9.0;
    cffs->imaginary_E        = (F2u*interpol_e_u_im*2.0 + F2d*interpol_e_d_im)/9.0;

    cffs->real_Etilde        = hA*(4.0*interpol_pion_pole_u_re + interpol_pion_pole_d_re)/9.0;
    cffs->imaginary_Etilde   = 0.0; //modelled directly after real-valued pion-pole
  }

  return;
}


void calcBH(struct _beam beam,struct _kinematics kinematics, 
	      struct _form_factors form_factors, struct _exclusive_vectors * vectors,
	    double E_p, double E_scatter, double E_gamma, double nu,
	    double qmod, double cos_t_gg, double tau, double lam1, double lam2, double produced_g_phi, double * dsBH){

  fprintf(stderr,"CalcBH: %f %f %f %f %f %f %f %f %f\n",
	  E_p,E_scatter,E_gamma,nu,qmod,cos_t_gg,tau,lam1,lam2);


  double k1pl = beam.energy*(Mp + E_p) - vectors->incidentLepton.Dot(vectors->outgoingTarget);
  double k2pl = E_scatter*(Mp + E_p) - vectors->outgoingLepton.Dot(vectors->outgoingTarget);
  double Mlep,Mlep2;

  // Choose your lepton!
  if(beam.particle==1) {Mlep = Mele; Mlep2 = Mele2;}
  if(beam.particle==2) {Mlep = Mmu; Mlep2 = Mmu2;}

#ifdef __TESTING_BH__

  fprintf(stderr,"Incident Lepton 3 vector: %lf %lf %lf\n",vectors->incidentLepton(0),vectors->incidentLepton(1),vectors->incidentLepton(2));
  fprintf(stderr,"Outgoing Target 3 vector: %lf %lf %lf\n",vectors->outgoingTarget(0),vectors->outgoingTarget(1),vectors->outgoingTarget(2));
  fprintf(stderr,"beam.energy,Mp,E_p,vectors->incidentLepton.Dot(vectors->outgoingTarget)");
  fprintf(stderr,"%lf %lf %lf %lf\n",beam.energy,Mp,E_p,vectors->incidentLepton.Dot(vectors->outgoingTarget));
  fprintf(stderr,"E_scatter,vectors->outgoingLepton.Dot(vectors->outgoingTarget)");
  fprintf(stderr,"%lf %lf\n",E_scatter,vectors->outgoingLepton.Dot(vectors->outgoingTarget));
#endif

  // Split these calcs up for re-use
  double summand1 = -lam1/lam2 - lam2/lam1;
  double summand2 = ( Mlep2*kinematics.del2*pow((1.0/lam1 - 1.0/lam2),2) ) / 2.0;
  double summand3 = kinematics.Q2*kinematics.del2/(2.0*lam1*lam2);
  double summand4 = pow( (k1pl/lam2 - k2pl/lam1)*Mlep/Mp ,2)/(2.0*(1.0-tau));
  double summand5 = ( tau/(1.0-tau) )*( pow(k1pl,2) + pow(k2pl,2))/(lam1*lam2);

  // similar to summand2, except for this sign change----|-|-------
  double summand6 = ( Mlep2*kinematics.del2*pow((1.0/lam1 + 1.0/lam2),2) ) / 2.0; 
  double summand7 = 2.0*Mlep2*(1.0/lam1 - 1.0/lam2)*(2.0 - Mlep2*(1.0/lam1 - 1.0/lam2));
  
  #ifdef __TESTING_BH__
  fprintf(stderr,"lam1,lam2,Mlep2,k1pl,k2pl,tau\n");
  fprintf(stderr,"%lf %lf %lg %lf %lf %lf\n",lam1,lam2,Mlep2,k1pl,k2pl,tau);
  #endif
  

  double Y1 = summand1 - summand2 + summand3 - summand4 - summand5;
  double Y2 = summand1 - summand6 + summand3 + summand4 + summand5 +summand7;
  
  double GE = form_factors.F1p + tau*form_factors.F2p; // tau is negative...
  double GM = form_factors.F1p + form_factors.F2p;

#ifdef __TESTING__

  fprintf(stderr,"Y1 %lf; Y2 %lf\n",Y1,Y2);
  
#endif


  double T2 = ((8.0*Mp2)/pow(kinematics.del2,2))*(pow(GE,2)*Y1 + tau*pow(GM,2)*Y2);

  
  T2 = T2*(E_scatter/beam.energy) / (16.0*Mp*pow(PI,2));

#ifdef __TESTING__
  
  fprintf(stderr,"GE %lf; GM %lf\n",GE,GM);
  fprintf(stderr,"T2 %lf\n",T2);
#endif

  
  *dsBH = T2*E_gamma/(Mp + nu - qmod*cos_t_gg);
 
  fprintf(stderr,"Pros: %lf %lf\n",k1pl,k2pl);
  fprintf(stderr,"dsBH from calcBH: %lf\n",*dsBH);

  //  double temp = calcBH1(beam,kinematics,form_factors,produced_g_phi,tau,&temp);

  //  fprintf(stderr,"dsBH from CalcBH1: %lg\n",temp);
  return;
}	     

void calcDVCS(int gpdmod, int target, struct _beam beam, 
	      struct _kinematics * kinematics, struct _form_factors form_factors, 
	      double del2min, double del2max, double lam1, double lam2, double tau, double phip,
	      double *dsDVCS, double *dsINT){

#ifdef __TESTING__

  fprintf(stderr,"Call to calcDVCS\n");
  fprintf(stderr,"gpd model %d; target %d; beam charge %d; beam helicity %f; beam particle %d\n",
	  gpdmod,target,beam.charge,beam.helicity,beam.particle);
  fprintf(stderr,"Kinematics: Q2 %f; xB %f; t %f; y %f\n",kinematics->Q2,kinematics->xB,kinematics->del2,kinematics->y);
  fprintf(stderr,"Form Factors: F1p %f; F2p %f; F1n %f; F2n %f\n",form_factors.F1p,form_factors.F2p,form_factors.F1n,form_factors.F2n);
  fprintf(stderr,"              GE_p %f; GM_p %f; GE_n %f GM_n %f\n",form_factors.GE_p,form_factors.GM_p,form_factors.GE_n,form_factors.GM_n);
  fprintf(stderr,"Other values: del2min %f;del2max %f; lam1 %f; lam2 %f; tau %f; phip %f\n",
	  del2min,del2max,lam1,lam2,tau,phip);
#endif
  // Get the relevant CFF values
  struct _compton_form_factors cffs;
  calcCFFs(gpdmod,target, *kinematics, form_factors,&cffs);
  
  // Some handy calculation shorthand
  double calc_sh_1 = 1 - del2min/kinematics->del2;
  double xB2 = kinematics->xB*kinematics->xB;
  double calc_sh_2 = pow(2 - kinematics->xB,2);
  double cy2 = 2.0 - 2.0*kinematics->y + pow(kinematics->y,2);
  double kfac = sqrt((-kinematics->del2)*(1 - kinematics->xB)*(1 - kinematics->y)*calc_sh_1/kinematics->Q2);


  // DVCS Contribution
  double a1 = pow(cffs.real_H,2) + pow(cffs.imaginary_H,2) + pow(cffs.real_Htilde,2) + pow(cffs.imaginary_Htilde,2);

  double a2 = 2 * (cffs.real_H*cffs.real_E + cffs.imaginary_H*cffs.imaginary_E +
		  cffs.real_Htilde*cffs.real_Etilde + cffs.imaginary_Htilde*cffs.imaginary_Etilde);

  double a3 = pow(cffs.real_E,2) + pow(cffs.imaginary_E,2);
  double a4 = pow(cffs.real_Etilde,2) + pow(cffs.imaginary_Etilde,2);

  // C-functions
  double C_DVCS = ( 4*(1 - kinematics->xB)*a1 - a2*xB2 - (xB2 + calc_sh_2*tau)*a3 - xB2*tau*a4 ) / calc_sh_2;
  double C_DVCS_eff = -kinematics->xB*C_DVCS;

  // Fourier coefficients
  double c0unp_DVCS = 2.0*cy2*C_DVCS;
  double c1unp_DVCS = 8.0*((2-kinematics->y)/(2-kinematics->xB))*C_DVCS_eff;

  // Scattering amplitude
  double T_DVCS = (c0unp_DVCS + kfac*c1unp_DVCS*cos(phip))/(pow(kinematics->y,2)*kinematics->Q2);

  double prefac = kinematics->xB*pow(kinematics->y,2)/
    (16*pow(PI,2)*pow(kinematics->Q2,2) * sqrt(1 + 4*Mp2*pow(kinematics->xB,2)/kinematics->Q2));
  
#ifdef __TESTING__

  fprintf(stderr,"C_DVCS %lf C_DVCS_eff %lf T_DVCS %lf prefac %lg\n",
	  C_DVCS,C_DVCS_eff,T_DVCS,prefac);
#endif

  // Cross section contribution!
  *dsDVCS = prefac*T_DVCS;
  
  // Interference contribution
  // Propagators
  double P1 = -2*lam1/kinematics->Q2;
  double P2 = 2*lam2/kinematics->Q2;

  // C functions
  double C_I_re = form_factors.F1p*cffs.real_H + 
    (kinematics->xB/(2-kinematics->xB))*(form_factors.F1p+form_factors.F2p)*cffs.real_Htilde -
    tau*form_factors.F2p*cffs.real_E;

  double C_I_im = form_factors.F1p*cffs.imaginary_H + 
    (kinematics->xB/(2-kinematics->xB))*(form_factors.F1p+form_factors.F2p)*cffs.imaginary_Htilde -
    tau*form_factors.F2p*cffs.imaginary_E;

  double C_I_re_eff = -kinematics->xB * C_I_re;
  double C_I_im_eff = -kinematics->xB* C_I_im;
  
  double RE2 = (kinematics->xB/(2-kinematics->xB))*(cffs.real_H + cffs.real_E) + cffs.real_Htilde;

  double b1 = (2-kinematics->xB)*(1-kinematics->y) - (1-kinematics->xB)*pow((2-kinematics->y),2)*calc_sh_1;
  double b2 = (1 - kinematics->y)*kinematics->xB*(form_factors.F1p + form_factors.F2p);

  double c0unp_I = -8*(2 - kinematics->y)*(b1*C_I_re - b2*RE2);
  double c1unp_I = -8*cy2*C_I_re;
  double c2unp_I = -16*((2-kinematics->y)/(2 - kinematics->xB))*C_I_re_eff;

  double s1unp_I = 8*kinematics->y*(2 - kinematics->y)*C_I_im;
  double s2unp_I = 16*(kinematics->y/(2 - kinematics->xB))*C_I_im_eff;

  double T_I = (kinematics->del2/kinematics->Q2)*c0unp_I + 
    kfac*(c1unp_I*cos(phip) + beam.helicity*s1unp_I*sin(phip)) + 
    pow(kfac,2)*(c2unp_I*cos(2*phip) + beam.helicity*s2unp_I*sin(2*phip));

#ifdef __TESTING__

  fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",RE2,b1,b2,c0unp_I,s1unp_I,s2unp_I);
  fprintf(stderr,"C_I_re %lf; T_I %lf\n",C_I_re,T_I);
#endif


  fprintf(stderr,"Lepton propogators from the interference term: %lf %lf\n",P1,P2);

  T_I = T_I/(kinematics->xB*pow(kinematics->y,3)*P1*P2*(-kinematics->del2));

  *dsINT = prefac*T_I;

  return ;
    
}
  

// Read in the pdf values from the look up table into memory
int readTab(int gpdmod,
	    double * re_h_u, double * re_h_d, 
	    double * re_del_u, double * re_del_d, 
	    double * re_e_u, double * re_e_d, 
	    double * re_pion_pole_u, double * re_pion_pole_d,
	    double * im_h_u, double * im_h_d, 
	    double * im_del_u, double * im_del_d,
	    double * im_e_u, double * im_e_d, 
	    double * im_pion_pole_u, double * im_pion_pole_d){

  FILE * lkup = fopen("gpd.dat","r");
  char buf[2048];
  
  int lineNo=0;
  int targetLine,entry;
  int ii;
  int flag = 0;
  int linecnt[51];

  // Populate the E-tilde pion pole info
  while(fgets(buf,sizeof(buf)-1,lkup) && lineNo <=50){
      
      sscanf(buf,"%*f %*f %*f %*f %lf %lf",re_pion_pole_u + lineNo, re_pion_pole_d + lineNo);
      lineNo++;
  }

  rewind(lkup);
  lineNo = 1;

  // Since the pion-pole is real valued, set imag. E-tilde to 0    
  for(ii = 0; ii < 51; ii++){
    im_pion_pole_u[ii] = 0.0;
    im_pion_pole_d[ii] = 0.0;
  }

  // Fill real cff components
  if (gpdmod == 1) targetLine = 1;
  else if (gpdmod == 2) targetLine = 52;
  else if (gpdmod == 3) targetLine = 103;
  else if (gpdmod == 4) targetLine = 154;
  else if (gpdmod == 5) targetLine = 205;
  else{
    fprintf(stderr,"Your gpd model was out of range: %d\n",gpdmod);
    exit(EXIT_FAILURE);
  }

  // Read in the real cff values from the file, taking only the values that are needed
  while(fgets(buf,sizeof(buf)-1,lkup)){
    
    #ifdef __1TESTING1__
    fprintf(stderr,"target %d lineNo %d flag %d entry %d\n",targetLine, lineNo, flag, entry);
    #endif
    if (lineNo++ != targetLine && flag == 0) continue;
    
    flag = 1;
    entry = lineNo - targetLine - 1;
    linecnt[entry] = lineNo;

    if(gpdmod == 1)
      sscanf(buf,"%lf %lf %lf %lf %*f %*f", re_h_u + entry,re_h_d + entry,re_del_u + entry,re_del_d + entry);
    else if(gpdmod == 2 || gpdmod == 3)
      sscanf(buf,"%lf %lf %lf %lf", re_h_u + entry,re_h_d + entry,re_del_u + entry,re_del_d + entry);
    else 
      sscanf(buf,"%lf %lf %lf %lf %lf %lf", 
	     re_h_u + entry,re_h_d + entry,re_del_u + entry,re_del_d + entry,re_e_u + entry,re_e_d + entry);

    if (entry ==  50) {flag = 0; break;}
  }

  if ( gpdmod <= 3){
    
    for(ii = 0; ii < 51; ii++){
      re_e_u[ii] = re_h_u[ii];
      re_e_d[ii] = re_h_d[ii];
    }
  }

  // Fill imaginary cff components
  if (gpdmod == 1) targetLine = 1+255;
  else if (gpdmod == 2) targetLine = 52+255;
  else if (gpdmod == 3) targetLine = 103+255;
  else if (gpdmod == 4) targetLine = 154+255;
  else if (gpdmod == 5) targetLine = 205+255;
  else{
    fprintf(stderr,"Your gpd model was out of range: %d\n",gpdmod);
    exit(EXIT_FAILURE);
  }

     
  while(fgets(buf,sizeof(buf)-1,lkup)){
    
    if (lineNo++ != targetLine && flag == 0) continue;

    flag = 1;
    entry = (lineNo  - targetLine) - 1;
    linecnt[entry] = lineNo;
    
    if(gpdmod == 1)
      sscanf(buf,"%lf %lf %lf %lf %*f %*f", im_h_u + entry,im_h_d + entry,im_del_u + entry,im_del_d + entry);
    else if(gpdmod == 2 || gpdmod == 3)
      sscanf(buf,"%lf %lf %lf %lf", im_h_u + entry,im_h_d + entry,im_del_u + entry,im_del_d + entry);
    else 
      sscanf(buf,"%lf %lf %lf %lf %lf %lf", 
	     im_h_u + entry,im_h_d + entry,im_del_u + entry,im_del_d + entry,im_e_u + entry,im_e_d + entry);
    
    if (entry == 50) {flag = 0; break;}
  }
  
  if ( gpdmod <= 3){

    for(ii = 0; ii < 51; ii++){
      im_e_u[ii] = im_h_u[ii];
      im_e_d[ii] = im_h_d[ii];

    }
  }

  fclose(lkup);

  // Sanity check
#ifdef __1TESTING1__
  fprintf(stderr,"pdf table\t gpd model %d\n",gpdmod);
  fprintf(stderr,"Real Values:\n");
  for(ii = 0; ii < 51; ii++)
    fprintf(stderr,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",linecnt[ii],
	    re_h_u[ii],re_h_d[ii],re_del_u[ii],re_del_d[ii],re_e_u[ii],re_e_d[ii],re_pion_pole_u[ii],re_pion_pole_d[ii]);
  fprintf(stderr,"\nImaginary Values:\n");
  for(ii = 0; ii < 51; ii++)
    fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf %lf\n",
	    im_h_u[ii],im_h_d[ii],im_del_u[ii],im_del_d[ii],im_e_u[ii],im_e_d[ii],im_pion_pole_u[ii],im_pion_pole_d[ii]);
  fprintf(stderr,"END\n");
#endif
  
  return 0;
}

		       
int bhdvcs(int xsecopt, int xsectyp, int gpdmod, int target, 
	   struct _beam beam,struct _kinematics *kinematics, struct _exclusive_vectors * vectors, 
	   double scattered_e_phi, double produced_g_phi,
	   double *dsBH, double *dsDVCS, double *dsINT){

  double Mlep,Mlep2;

  if(beam.particle==1){Mlep = Mele; Mlep2 = Mele2;}
  if(beam.particle==2){Mlep = Mmu; Mlep2 = Mmu2;}

  double Mp2 = Mp*Mp;

  double xmin = kinematics->Q2/(2*Mp*beam.energy);
  double xmax = 1.0;

#ifdef __COMPASS__
  
  produced_g_phi = produced_g_phi - PI;

#endif

#ifdef __1TESTING__

  fprintf(stderr,"GPD Model is %d; XS is type %d; target is %d\n",gpdmod,xsectyp,gpdmod);
  fprintf(stderr,"beam profile: particle %d; helicity %d; charge %d\n",beam.particle,beam.helicity,beam.charge);
  fprintf(stderr,"Kinematics: Q2 %f; xB %f; t %f\n",kinematics->Q2,kinematics->xB,kinematics->del2);
  fprintf(stderr,"e_phi %f; g_phi %f\n",scattered_e_phi,produced_g_phi);
#endif

  struct _form_factors form_factors;
  //   Sanity check
  if( kinematics->xB > xmax || kinematics->xB < xmin){
    fprintf(stderr,"xMax: %lf\txMin: %lf  ---  Bjorken X variable is out of range: %lf\n", xmax,xmin,kinematics->xB);
    return -1;
  }

  // Calculate other kin vars
  double nu = kinematics->Q2/(2*Mp*kinematics->xB);
  double yB = nu / beam.energy;
  kinematics->y = yB;
  if(yB < 0 || yB > 1){
    fprintf(stderr,"Bjorken y out of range: %lf\n",yB);
    return -1;
  }
  
  double E_scatter = beam.energy - nu;
  if(E_scatter < 0){
    fprintf(stderr,"Scattered particle energy was out of range: %lf\n",E_scatter);
    return -1;
  }
  
  double W2 = Mp2 + 2*Mp*nu - kinematics->Q2;
  if (W2 < 0){
    fprintf(stderr,"W2 out of range: %lf\n",W2);
    return -1;
  }

  // Convenience terms
  double W = sqrt(W2);
  double qmod = sqrt(nu*nu + kinematics->Q2);


  double E1com = Mp*(Mp+nu)/W;
  double P1com = Mp*qmod/W;
  
  double E2com = (W2 + Mp2)/(2.0*W);
  double P2com = (W2 - Mp2)/(2.0*W);
  
  // Min and Max mom transfer
  double del2min = 2.0*(Mp2 - E1com*E2com + P1com*P2com);
  double del2max = 2.0*(Mp2 - E1com*E2com - P1com*P2com);
  fprintf(stderr,"Del2min is %lf\n",del2min);
  if(kinematics->del2 > del2min || kinematics->del2 < del2max){
    fprintf(stderr,"Del2 is out of range: %lf\n",kinematics->del2);
    return -1;
  }

  // cos(theta) scattered lepton
  double costel = 1.0 - (kinematics->Q2/(2.0*beam.energy*E_scatter)); 
  if(costel > 1.0){
    fprintf(stderr,"costel out of range: %lf\n",costel);
    return -1;
  }
  
  // sin theta scattered lepton
  double sintel = sqrt(1 - costel*costel); 

  
  // Energy of the scattered proton
  double E_p = Mp - kinematics->del2/(2.0*Mp);
  
  // Energy of the produced photon 
  double E_gamma = nu + kinematics->del2/(2.0*Mp);
  
  // Convenience term
  double tau = kinematics->del2/(4.0*Mp2);

  // Fill the three-vector directions of the various particles
  vectors->incidentLepton.SetX(0.0);
  vectors->incidentLepton.SetY(0.0);
  vectors->incidentLepton.SetZ(beam.energy);
  
  vectors->outgoingLepton.SetX(E_scatter*sintel);
  vectors->outgoingLepton.SetY(0.0);
  vectors->outgoingLepton.SetZ(E_scatter*costel);

  #ifdef __TESTING__
  
  fprintf(stderr,"Outgoing Lepton: %lf %lf %lf\n",vectors->outgoingLepton(0),vectors->outgoingLepton(1),vectors->outgoingLepton(2));
  #endif

  vectors->virtualPhoton = vectors->incidentLepton - vectors->outgoingLepton;

  // cos theta virtual photon
  double cos_t_Vq = vectors->virtualPhoton(2)/qmod;
  // sin theta virtual photon
  double sin_t_Vq = sqrt(1.0-(cos_t_Vq*cos_t_Vq));
  
  // cos theta gamma gamma*
  double cos_t_gg = (2.0*E_gamma*(Mp+nu) + kinematics->Q2 - 2.0*Mp*nu)/(2.0*E_gamma*qmod);
  // sin theta gamma gamma*
  double sin_t_gg = sqrt(1.0 - (cos_t_gg*cos_t_gg));
  
  // Convenience terms
  double facX = E_gamma*sin_t_gg*cos(produced_g_phi);
  double facY = E_gamma*sin_t_gg*sin(produced_g_phi);
  double facZ = E_gamma*cos_t_gg;

  // Fill outgoing photon vector
  vectors->outgoingPhoton.SetX(facX*cos_t_Vq - facZ*sin_t_Vq);
  vectors->outgoingPhoton.SetY(facY);
  vectors->outgoingPhoton.SetZ(facX*sin_t_Vq + facZ*cos_t_Vq);

#ifdef __TESTING__
  fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf\n",facX,facY,facZ,E_gamma,cos_t_Vq,sin_t_Vq,produced_g_phi);
  fprintf(stderr,"Virtual Photon 3 vector: %lf %lf %lf\n",vectors->virtualPhoton(0),vectors->virtualPhoton(1),vectors->virtualPhoton(2));
  fprintf(stderr,"Outgoing Photon 3 vector: %lf %lf %lf\n",vectors->outgoingPhoton(0),vectors->outgoingPhoton(1),vectors->outgoingPhoton(2));
#endif
  // Calcualte outgoing proton vector
  vectors->outgoingTarget = vectors->virtualPhoton - vectors->outgoingPhoton;

  // cos theta incident lepton outgoing photon
  double cos_t_eg = ( vectors->incidentLepton.Dot(vectors->outgoingPhoton) )/ (beam.energy*E_gamma);
  double te_t_eg = acos(cos_t_eg);
  
  // cos theta scattered lepton outgoing photon
  double cos_t_scat_g = 
    ( vectors->outgoingLepton.Dot(vectors->outgoingPhoton) )/ (E_scatter * E_gamma);
  double te_t_scat_g = acos(cos_t_scat_g);

  double lam1 = E_gamma*( 2*beam.energy*pow(sin(te_t_eg/2.0),2)   + Mlep2*cos_t_eg/(2.0*beam.energy));
  double lam2 = E_gamma*( 2*E_scatter*sq(sin(te_t_scat_g/2.0)) + Mlep2*cos_t_scat_g/(2.0*E_scatter));

  #ifdef __TESTING__
  fprintf(stderr,"E_gamma %lf; te_t_eg %lf; te_scat_e_g %lf; E_scatter %lf\n",E_gamma,te_t_eg,te_t_scat_g,E_scatter);
  #endif

  double scattered_p_phi = -produced_g_phi + PI;

  calcFFs(kinematics->del2,&form_factors);
  
  if( xsecopt & 1 )  calcBH( beam, *kinematics, form_factors, vectors,
			     E_p, E_scatter, E_gamma, nu,
			     qmod, cos_t_gg, tau, lam1, lam2, produced_g_phi,dsBH);

  if( xsecopt & 2 )  calcDVCS( gpdmod, target, beam, kinematics, form_factors, 
			       del2min, del2max, lam1, lam2, tau, scattered_p_phi, dsDVCS, dsINT);

  double cjakob = (pow(kinematics->y,2)/(1-kinematics->y))*(W2 - Mp2)/(4.0*kinematics->Q2*qmod*pow(E_gamma,2));

  if(xsectyp==1)
    *dsBH = (*dsBH)*cjakob;
  else{
    *dsDVCS = *dsDVCS/cjakob;
    *dsINT = *dsINT/cjakob;
  }

  // Convert into inverse picobarn uits - coeff is 1 pb^-1
  *dsBH = coeff*(*dsBH);
  *dsDVCS = coeff*(*dsDVCS);
  *dsINT = -beam.charge* coeff*(*dsINT);

  return 0;
}


#ifdef __TESTING1__
int main(){

  struct _beam beam;
  struct _kinematics kinematics;
  struct _exclusive_vectors vectors;
  
  double scattered_e_phi = 0;
  double produced_g_phi = 0.;
  
  double dsBH, dsDVCS,dsINT;

  // Test kinematics
  beam.charge = 1;
  beam.helicity = 0;
  beam.particle = 1;
  beam.energy = 27.56;

  kinematics.Q2 = 2.4;
  kinematics.xB = 0.1;
  kinematics.del2 = -0.1;

  bhdvcs(3,1,1,1,beam,&kinematics,&vectors,
	 scattered_e_phi,produced_g_phi,
	 &dsBH,&dsDVCS,&dsINT);

  fprintf(stdout,"BH Cross Section was %lf\n",dsBH);
  fprintf(stdout,"DVCS Cross Section was %lf\n",dsDVCS);
  fprintf(stdout,"Interference Cross Section was %lf\n",dsINT);

  return 0;
}
#endif 
