#include "bhdvcs.h"

double calcBH1(struct _beam beam, struct _kinematics ks, 
	       struct _form_factors ff, double phig, double tau,double *dsBH){

  double eps = 2*ks.xB*Mp/sqrt(ks.Q2);
  double top = -ks.Q2*(2*(1-ks.xB)*(1-sqrt(1+sq(eps)))+sq(eps));
  double bottom = 4*ks.xB*(1-ks.xB)+sq(eps);
 
  fprintf(stderr,"eps: %lf\n",eps);
  fprintf(stderr,"top: %lf\n",top);
  fprintf(stderr,"bottom: %lf\n",bottom);

  double delmin = -1*(ks.Q2*2*(1-ks.xB)*(1-sqrt(1+sq(eps)))+sq(eps))
                                 /(4*ks.xB*(1-ks.xB)+sq(eps));
 
  double delmint =  -(Mp2*sq(ks.xB)/(1-ks.xB+ks.xB*Mp2/ks.Q2));

  fprintf(stderr,"Delmin test %g %lf\n",delmin,delmint);

  double Kfac = sqrt(-(ks.t/ks.Q2)*(1-ks.xB)*(1 - ks.y - sq(ks.y)*sq(eps)/4)*(1-delmint/ks.t)*(sqrt(1+sq(eps)) 
											       + (4*ks.xB*(1-ks.xB)+sq(eps))/(4*(1-ks.xB))*(ks.t-delmint)/ks.Q2)); //nan

  double fac1 = -(ks.del2/ks.Q2);
  double fac2 = 1-ks.xB;
  double fac3 = 1-ks.y-sq(ks.y)*sq(eps)/4;
  double fac4 = 1 - delmin/ks.t;
  double fac5 = sqrt(1+sq(eps)) + (4*ks.xB*(1-ks.xB)+sq(eps))*(ks.del2-delmin)/(4*(1-ks.xB)*ks.Q2);

  double Kfac2 = fac1*fac2*fac3*fac4*fac5;
  
//   fprintf(stderr,"*****************************************************\n");
//   fprintf(stderr,"*****************************************************\n");
//   fprintf(stderr,"%lf %lf %lf %lf %lf\n",fac1,fac2,fac3,fac4,fac5);
//   fprintf(stderr,"*****************************************************\n");
//   fprintf(stderr,"*****************************************************\n");
  Kfac = sqrt(Kfac2);

  double kdel = -ks.Q2/(2*ks.y*(1+eps)) * 
    (1 + 2*Kfac*cos(phig) - (ks.t/ks.Q2)*(1-ks.xB*(2-ks.y)+ks.y*pow(eps,2)/2) + ks.y*pow(eps,2)/2); //nan

  double prop1 = (ks.Q2+2*kdel)/ks.Q2; 
  double prop2 = (-2*kdel+ks.t)/ks.Q2;
  
  fprintf(stderr,"Lepton propogators from calcBH1: %lf %lf\n",prop1,prop2);


  double sqffdiff = sq(ff.F1p) + tau*sq(ff.F2p); //GE^2
  double sqffadd = sq(ff.F1p+ff.F2p); //GM^2

//   fprintf(stderr,"*****************************************************\n");
//   fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf\n",
// 	  Kfac2,Kfac,kdel, prop1,prop2,sqffdiff,sqffadd);
//   fprintf(stderr,"*****************************************************\n");


  double summand1 = 8*sq(Kfac)*((2+3*sq(eps))*(ks.Q2/ks.t)*sqffdiff + 2*sq(ks.xB)*sqffadd) ;

  double summand2 = sq(2-ks.y)*((2+sq(eps))*((sq(ks.xB)/tau)*sq(1 + ks.t/ks.Q2) 
					     + 4*(1-ks.xB)*(1+ks.xB*ks.t/ks.Q2))*sqffdiff) ;

  double summand3 = 4*sq(ks.xB)*(ks.xB + 
				 (1-ks.xB+sq(eps)/2)*sq(1-ks.t/ks.Q2) - ks.xB*(1-2*ks.xB)*sq(ks.t/ks.Q2)*sqffadd);

  double summand4 = 8*(1+sq(eps))*(1-ks.y-sq(eps*ks.y)/4)*(2*sq(eps)*(1-tau)*sqffdiff - sq(ks.xB)*sq(1-ks.t/ks.Q2)*sqffadd);

  double c0_BH = summand1+summand2+summand3+summand4;
    
  fprintf(stderr,"BH C0: %lf %lf %lf %lf %lf\n",summand1,summand2,summand3,summand4,c0_BH);
		 

  fac1 =     8*Kfac*(2-ks.y);

  fac2 = ((4*sq(ks.xB)*Mp2/ks.t - 2*ks.xB - sq(eps))*sqffdiff + 2*sq(ks.xB)*(1-(1-2*ks.xB)*ks.t/ks.Q2*sqffadd));
				  
  double c1_BH = fac1*fac2;
  
  fprintf(stderr,"BH C1: %lf %lf %lf\n",fac1,fac2,c1_BH);
  
  double c2_BH = 8*sq(ks.xB)*Kfac2*(((1/tau)*sqffdiff)+2*(sqffadd));
  
  // Returns eqn. 25 from BMK
  fprintf(stderr,"calcBH1: %lf %lf %lf\n",c0_BH,c1_BH,c2_BH);
  //  fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",eps,delmint,Kfac2,Kfac,kdel,prop1,prop2,sqffdiff,sqffadd);
  fprintf(stderr,"Cross section factors:\n");
  fprintf(stderr,"%lf %lf %lf %lf\n",sq(ks.xB*ks.y),1+sq(eps),ks.t,prop1*prop2);
  *dsBH = 1/(sq(ks.xB*ks.y)*sq(1+sq(eps))*ks.t*prop1*prop2) * (c0_BH + c1_BH*cos(phig) + c2_BH*cos(phig));
  fprintf(stderr,"%lf\n",*dsBH);
  fprintf(stderr,"%lg\n",pow(ALPHA_EM,3)*ks.xB*ks.y/(16*sq(PI)*ks.Q2*sqrt(1+sq(eps))));
  return *dsBH = pow(ALPHA_EM,3)*ks.xB*ks.y/(16*sq(PI)*ks.Q2*sqrt(1+sq(eps)))*(*dsBH);
  
}
