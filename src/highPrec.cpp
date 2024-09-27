#include "../include/highPrec.h"

using namespace std;

PrecFloat ExpEiComplexSum(PrecFloat MOD, PrecFloat PH, PrecFloat s,  bool MODE) {

  PrecFloat gamma= precEuler();
  PH= atan(PH);
  PrecFloat RE_Log= log(MOD);
  PrecFloat IM_Log= PH;


  
  //evaluate sum_k=1,infty  z^k /( k * k! )

  PrecFloat START_SUM, SUM_OLD;
  if(MODE==0) START_SUM= (RE_Log+gamma)*sin(s) -IM_Log*cos(s);
  else START_SUM= (RE_Log+gamma)*cos(s) +IM_Log*sin(s);

  PrecFloat SUM=START_SUM;
  SUM_OLD=0.0;

  
  bool Arg_zero = (MOD==0);
  
  if(!Arg_zero) {
  
    bool converged=false;

    unsigned long int k=1;

    PrecFloat MOD_TERM_k= MOD;
    SUM+=MOD*((MODE==0)?sin(s- PH):cos(s - PH));
    k++;
   
    while( ! converged) {

      if(k%30==0) SUM_OLD= SUM;

      MOD_TERM_k *=  MOD*(k-1)/(PrecFloat(k)*PrecFloat(k));
      SUM += MOD_TERM_k*((MODE==0)?sin(s- k*PH):cos(s - k*PH));
    
  
      //cout<<"SUM: k="<<k<<" : "<< SUM<<endl;

      k++;
      if(SUM_OLD == SUM) converged=true;
      
    }
  }

  return SUM;
  
}

PrecFloat Erfi( PrecFloat x) {

  PrecFloat preF= 2/sqrt(precPi());
  bool CONVERGED=false;
  PrecFloat FACT=x;
  PrecFloat SUMM= FACT;
  PrecFloat SUMM_OLD=SUMM;
  unsigned long int n=1;
  while(!CONVERGED) {
    FACT *= x*x/n;
    SUMM += FACT/(2*n+1);
    if(SUMM==SUMM_OLD) CONVERGED=true;
    SUMM_OLD=SUMM;
    n++;
    
  }

  //cout<<"nsteps to evaluate Erfi(x="<<x<<"): "<<n<<endl<<flush;
  return preF*SUMM;
}



PrecFloat DawsonF(PrecFloat x) {

  int sign= ((x >= 0)?1:-1);
  x= abs(x);
  
  bool MODE=0;
  bool CONVERGED=false;
  PrecFloat SUMM, SUMM_OLD, FACT;
  unsigned long int n=1;
   
  if( (x<=18.93)) { //use power expansion around x=0

    SUMM= Erfi(x)*exp(-x*x)*sqrt( precPi())/2;
    /*
    FACT= x;
    SUMM=FACT;
    SUMM_OLD=SUMM;
    while(!CONVERGED) { // sum (-1)^n * 2^n * x^(2n+1) /(2n+1)!!
      
      FACT *= -2*x*x*(2*n+1);

      SUMM += FACT;
      if(SUMM==SUMM_OLD) CONVERGED=true;
      SUMM_OLD=SUMM;
      cout<<"MODE 0: nstep: "<<n<<" SUMM: "<<SUMM<<endl<<flush;
      n++;
      
    }
    */
    
  }
  else { //use power expansion around x=Infinity
    MODE=1;
    FACT= 1/(2*x);
    SUMM=FACT;
    SUMM_OLD=SUMM;
    while(!CONVERGED) { // \sum_n (2n-1)!! /( 2^(n+1) * x^(2n+1) )
      
      FACT /= 2*x*x/(2*n-1);
      SUMM+=FACT;
      if(SUMM==SUMM_OLD) CONVERGED=true;
      SUMM_OLD=SUMM;
      n++;
    }
    //cout<<"DawsonF(x="<<x<<"): "<<SUMM<<" nsteps: "<<n<<" MODE: "<<MODE<<endl<<flush;
  }
  
  return sign*SUMM;
 
}
