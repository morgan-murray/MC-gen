#include "bhdvcs.h"



int main(){


  struct _beam bm;
  struct _kinematics k;
  struct _exclusive_vectors v;

  double e_phi, g_phi;
  double dsDVCS, dsBH, dsINT;


  char fname[]="compass_kinrange_posmu_neghel.dat";
  FILE *f = fopen(fname,"w");

  bm.charge = 1;
  bm.helicity = -0.8;
  bm.particle = 1;
  bm.energy = 27.56;
  
  k.Q2=2.4;
  k.del2 = -0.1;
  k.xB = 0.05;
  k.t = k.del2;
  
  e_phi=0;

  double x;
  double t;
  double q;

  fprintf(f,"Approx. COMPASS Kinematic Range: 1.0 < Q2 < 8.0\n");
  fprintf(f,"                                 0.01 < -t < 0.2\n");
  fprintf(f,"                                 0.05 < xB < 0.15\n");
  fprintf(f,"Using mu+ beam, polarised to -0.8 helicity\n\n\n");
  fprintf(f,"Q2 t xB phi dsBH dsDVCS dsINT\n\n");
  
  for(q=1.0;q<8.0;q+=1){

    k.Q2=q;
    
    for(t = -0.01; t >= -0.2; t = t-0.01){
      
      k.del2=t;
      k.t=k.del2;
      
      
      for(x=0.05;x<=0.15;x+=0.01){
	
	k.xB=x;
	
	for(g_phi=-3.141;g_phi<=3.141;g_phi+=PI/30){
	  //	g_phi=0;
	  bhdvcs(3,1,4,1,bm,&k,&v,e_phi,g_phi,&dsBH,&dsDVCS,&dsINT);
	  
	  fprintf(f,"%lf %lf %lf %lf %lf %lf %lf\n",k.Q2,k.del2,k.xB,g_phi,dsBH,dsDVCS,dsINT);
	  
	  
	}    
      }
    }
  }
  fclose(f);
  
  return 0;
}
