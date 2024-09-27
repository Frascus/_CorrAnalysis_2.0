#include "../include/numerics.h"

using namespace std;


void D(int k) { 
  cout<<"D("<<k<<")"<<endl;
}

double eps(int k) {
  double epsilon= 1;
  for(int i=0;i<k;i++) epsilon/=10;

  return epsilon;
}
long long int fact(int n) {
  if(n==0) return 1;
  else return n*fact(n-1);
}


long long int ipow(int a,int n) {
  long long int result= a;
  if(a < 0 || n < 0) crash("ERR0R: ipow called with negative arguments");
  if(n==0) return 1;
  else return result*ipow(a,n-1);
}

double fpow(double a, int n) {
  double result = a;
  if(n<0) crash("ERROR: fpow called with negative arguments");
  if(n==0) return 1;
  else return result*fpow(a,n-1);
}

long long int BinomialCoeff(int n,int k) {

  vector<long long int> C(k+1,0);
   
    C[0] = 1;  // nC0 is 1 
  
    for (int i = 1; i <= n; i++) 
    { 
        // Compute next row of pascal triangle using 
        // the previous row 
        for (int j = min(i, k); j > 0; j--) 
            C[j] = C[j] + C[j-1]; 
    } 
    return C[k];   
}


double RatioPol(Vfloat &x_num, Vfloat &x_den, Vfloat &Num, Vfloat &Den, Vint &NumPow, Vint &DenPow, string MODE_POL) {

  if(x_num.size() == 0 && x_den.size() == 0) crash("RatioPol called with zero NUM and DEN");

  double NUM_VAL=0;
  double DEN_VAL=0;
  double VAL=0;

  if(MODE_POL=="SINGLE_RATIO") {

  if(x_num.size()==0) NUM_VAL=1.0;
  else {
    for(int i=0;i<(signed)x_num.size();i++) NUM_VAL += pow(x_num[i],NumPow[i])*Num[i];
  }
  if(x_den.size() ==0) DEN_VAL=1.0;
  else {

    for(int j=0; j<(signed)x_den.size();j++) DEN_VAL += pow(x_den[j],DenPow[j])*Den[j];
  }
  VAL= NUM_VAL/DEN_VAL;

  }
  else if(MODE_POL=="MULTIPLE_RATIO") {
    if((signed)x_num.size() != (signed)x_den.size()) crash("RatioPol called in MULTIPLE_RATIO MODE but x_num.size and y_num.size do not coincide");
    for(int i=0; i<(signed)x_num.size();i++) VAL+= pow(x_num[i],NumPow[i])*Num[i]/(pow(x_den[i],DenPow[i])*Den[i]); 
  }

  else crash("MODE: "+MODE_POL+" in RatioPol not yet implemented");


  return VAL;
}

void crash(string Message) {
  cout << Message <<endl;
  exit(-1);
  return ;
}


double F(double x, double a, int alpha) {
  
  if(alpha==0) return a-1.0/cosh((alpha+1)*x);
  else if(alpha==-1) return a -cosh(alpha*x);
  
  double return_value= a -cosh(alpha*x)/cosh((alpha+1)*x);
  return return_value;
}


void derivative(Vfloat &RES, Vfloat &INPUT, string MODE) {

  int size= INPUT.size();

  for(int t=0; t<size;t++) {
    int t1,t2;
    if(MODE=="BACKWARD") {t2=t; t1= (t2-1)%size;}
    if(MODE=="FORWARD") {t1=t; t2= (t1+1)%size;}
    if(MODE=="SYMMETRIC") {t1= (t-1)%size; t2=(t+1)%size;}

    double temp_res= INPUT[t2]-INPUT[t1];
    if(MODE=="SYMMETRIC") temp_res /= 2.0;
    RES.push_back(temp_res);


  }

  return;
}


double FTAN(double x1,  int t, int NT) {

  double dNT = (double)NT/2.0;
  double dt = (double)t;
  
  double RES=  +1.0*(dNT -dt +1.0)*tanh(x1*(dNT-dt+1.0)) -1.0*(dNT -dt)*tanh(x1*(dNT-dt));

  // cout<<"FTAN: "<<RES<<endl<<flush;

  return RES;

}

double FTAN_SYMM(double x1,  int t, int NT) {

  double dNT = (double)NT/2.0;
  double dt = (double)t;
  
  double RES=  +1.0*(dNT -dt +1.0)*tanh(x1*(dNT-dt+1.0)) -1.0*(dNT -dt-1.0)*tanh(x1*(dNT-dt-1.0));

  // cout<<"FTAN: "<<RES<<endl<<flush;

  return RES;


}



double Root_Brent(double R, int nt, int NT) {

  
   
  //cout<<"R: "<<R<<" t:"<<nt<<" T: "<<NT<<endl;
  if( R < 1 && nt < NT/2) return log(fabs(R)) ;
  if( R > 1 && nt >= NT/2) return log(1.0/R);
  
  if(NT%2 != 0) crash("Temporal lattice extent is not divisible by two!");
  NT= NT/2; //it is what enter the expression of the effective mass
  int alpha = nt-NT;
  if(alpha==0 && R<= 1) return acosh(1.0/R);
  double Precision = 1e-9;
  double delta = 0.001;
  //solve the equation cosh(meff(alpha+1))*corr_ratio = cosh(meff*alpha) using the Brent method!
  //initialize iteration
  double b=0;
  double a;
  a= min(fabs(1./(alpha+1)),fabs(1./alpha));

  while (F(a,R,alpha)*F(b,R,alpha) > 0) {a  +=  min(fabs(1./(alpha+1)),fabs(1./alpha));}
 
 



  if(fabs(F(a,R,alpha)) < fabs(F(b,R,alpha))) {double atemp=a; a=b; b=atemp;}

  double c=a;
  bool FLAG = true;
  double s=b;
  double d=0;

  
  
  while(F(s,R,alpha) !=0 && fabs(b-a)>= Precision*fabs((double)(b+a)/2.0) ) {

  
    if((F(a,R,alpha) != F(c,R,alpha)) && (F(b,R,alpha) != F(c,R,alpha))) {//inverse quadratic interpolation
      s= a*F(b,R,alpha)*F(c,R,alpha)/((F(a,R,alpha)-F(b,R,alpha))*(F(a,R,alpha)-F(c,R,alpha))) + b*F(a,R,alpha)*F(c,R,alpha)/((F(b,R,alpha)-F(a,R,alpha))*(F(b,R,alpha)-F(c,R,alpha))) + c*F(a,R,alpha)*F(b,R,alpha)/((F(c,R,alpha)-F(a,R,alpha))*(F(c,R,alpha)-F(b,R,alpha)));
    }
    else s= b-(F(b,R,alpha)*(b-a)/(F(b,R,alpha)-F(a,R,alpha)));

    double s1= (double)(3*a+b/4.0);
    if( (s < s1 || s> b) || (FLAG==true && fabs(s-b) >= (double)fabs(b-c)/2) || (FLAG==false && fabs(s-b) >= (double)fabs(c-d)/2) || (FLAG==true && fabs(b-c) < delta) || (FLAG==false && fabs(c-d) < delta)) {
      
      FLAG=true;
      s= (a+b)/2.0;
      
    }
    
    else FLAG= false;

    d= c;
    c= b;
    if (F(a,R,alpha)*F(s,R,alpha)<0) b=s;
    else a=s;

    if(fabs(F(a,R,alpha)) < fabs(F(b,R,alpha))) {double atemp=a; a=b; b=atemp;}

   
  }

  
  return s;
}



double Root_Brent_sinh(double R, int nt, int NT) {


 
  
   
  if(NT%2 != 0) crash("Temporal lattice extent is not divisible by two!");
  NT= NT/2; //it is what enter the expression of the effective mass
  int alpha = nt-NT;
  if(alpha==-1) return 0.0;

  auto F=[&](double x) {return R- sinh(x*alpha)/sinh(x*(alpha+1)); };
  double Precision = 1e-9;
  double delta = 0.001;
  //solve the equation sinh(meff(alpha+1))*corr_ratio = sinh(meff*alpha) using the Brent method!
  //initialize iteration
  double b=1e-8;
  double a;
  a= min(fabs(1./(alpha+1)),fabs(1./alpha));

  while (F(a)*F(b) > 0) {a  +=  min(fabs(1./(alpha+1)),fabs(1./alpha));}
 
 



  if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

  double c=a;
  bool FLAG = true;
  double s=b;
  double d=0;

  
  
  while(F(s) !=0 && fabs(b-a)>= Precision*fabs((double)(b+a)/2.0) ) {

  
    if((F(a) != F(c)) && (F(b) != F(c))) {//inverse quadratic interpolation
      s= a*F(b)*F(c)/((F(a)-F(b))*(F(a)-F(c))) + b*F(a)*F(c)/((F(b)-F(a))*(F(b)-F(c))) + c*F(a)*F(b)/((F(c)-F(a))*(F(c)-F(b)));
    }
    else s= b-(F(b)*(b-a)/(F(b)-F(a)));

    double s1= (double)(3*a+b/4.0);
    if( (s < s1 || s> b) || (FLAG==true && fabs(s-b) >= (double)fabs(b-c)/2) || (FLAG==false && fabs(s-b) >= (double)fabs(c-d)/2) || (FLAG==true && fabs(b-c) < delta) || (FLAG==false && fabs(c-d) < delta)) {
      
      FLAG=true;
      s= (a+b)/2.0;
      
    }
    
    else FLAG= false;

    d= c;
    c= b;
    if (F(a)*F(s)<0) b=s;
    else a=s;

    if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

   
  }

  
  return s;
}


double DoConstantFit(Vfloat &data, Vfloat &err) {

  double result=0;
  double weight=0;


  

  if (data.size() != err.size()) {
    cout<<"data size and err size do not match in double DoConstantFit(Vfloat&, Vfloat&)"<<endl;;
    cout<<"Printing info..."<<endl;
    printV(data, "data", 1);
    printV(err, "err", 1);
    crash("Abort.");
  }
  
  for(unsigned int k=0; k<data.size(); k++) {result += pow(1.0/err[k],2)*data[k]; weight += pow(1.0/err[k],2);}

  if(isnan(result) || isnan(weight)) { //fit without errors
    //cout<<"Switching from chi^2 minimization to least square minimization......"<<endl;
    result=0;
    weight=0;
    for(auto &meas: data) {result += meas; weight++;}
  }

  //normalize
  result /= weight;

  return result;
}

double lin_interpolator(double y1, double y2, double Dx1, double Dx2,double Dx) {

  double n1= y1*Dx2;
  double n2= y2*Dx1;

  double a0= (n1-n2)/(Dx2-Dx1);
  double a1 = (y1-a0)/Dx1;

  return a0+a1*Dx;
}


double quad_interpolator(double y1, double y2, double y3, double Dx1, double Dx2, double Dx3, double Dx) {


  
	double n1 = y1*Dx2*Dx3/( (Dx1-Dx2)*(Dx1-Dx3));
	double n2 = -y2*Dx1*Dx3/( (Dx1-Dx2)*(Dx2-Dx3));
	double n3 = y3*Dx1*Dx2/( (Dx1-Dx3)*(Dx2-Dx3));

	double prod= (Dx1-Dx2)*(Dx1-Dx3)*(Dx2-Dx3);

	double a0 = n1+n2+n3;
	double a1 = (pow(Dx3,2)*(y1-y2) + pow(Dx1,2)*(y2-y3) + pow(Dx2,2)*(y3-y1))/prod ;
	double a2 = (Dx3*(y2-y1) + Dx2*(y1-y3) + Dx1*(y3-y2))/prod;

	return a0 + a1*Dx + a2*pow(Dx,2);
  



}


void Print_To_File(const vector<string>& row_id, const VVfloat &data, string Path, string MODE, string Header) {

  ofstream Print;
  if(MODE=="APP") Print.open(Path, ofstream::app);
  else Print.open(Path, ofstream::out);

  Print.precision(16);

  if(data.size() == 0) crash("In Print_To_File an empty vectors has been provided");

  unsigned int rows= data[0].size();

  if((signed)row_id.size() != 0 && row_id.size() != rows) crash(" In Print_To_File data size and data_id size do not match"); 

  Print<<Header<<endl;
  for(unsigned int j=0; j < rows; j++) {

    Print<<j;
    if(row_id.size() != 0) Print<<setw(20)<<row_id[j];
    
    for(unsigned int i=0; i<data.size(); i++) {
      if(data[i].size() != rows) crash("In Print_To_File, number of measurements is non-constant over columns. Exiting...");
      Print<<setw(30)<<data[i][j];
    }
    Print<<endl;
  }
  Print.close();
  return;
}





void cascade_resize( vector<vector<double>>& arr, const Vint& A ) {

  if((signed)A.size() != 2) crash("In cascade_resize cannot resize vv<double> with a Vint of size != 2");
  arr.resize(A[0]);
  for( auto &a: arr) a.resize(A[1]);

  return;
  
}

void cascade_resize( vector<vector<vector<double>>>& arr,const Vint &A) {
 
  if((signed)A.size() != 3) crash("In cascade_resize cannot resize vvv<double> with a Vint of size != 3");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}

void cascade_resize( vector<vector<vector<vector<double>>>>& arr,const Vint &A) {
  
  if((signed)A.size() != 4) crash("In cascade_resize cannot resize vvvv<double> with a Vint of size != 4");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}


void cascade_resize( vector<vector<vector<vector<vector<double>>>>>& arr,const Vint &A) {

  if((signed)A.size() != 5) crash("In cascade_resize cannot resize vvvvv<double> with a Vint of size != 5");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}

bool Is_perfect_square(int x) {

  int sqrt_x = sqrt(x);
  if (sqrt_x*sqrt_x == x) return true;
  else return false;
}


int degeneracy(int m) {

  if(m==0) return 1;
  int deg_val=0;

  //return number of lattice sites with distance m>0 from origin.
  assert(m > 0);

  int x2_max= m;
  int x2_min;
  if(m%3==0) x2_min=m/3;
  else x2_min= floor( m/3 + 1) ;
  

  

  for(int x2=x2_min;x2<=x2_max;x2++) {
    if(Is_perfect_square(x2)) {  //x2 is a perfect square
      
      int y2_max= m-x2;
      int y2_min;
      if( (m-x2)%2 == 0) y2_min= (m-x2)/2;
      else y2_min = floor( (m-x2)/2 +1);
      
      
      for(int y2=y2_min;y2<=y2_max;y2++) {
	
	if(y2<=x2 && Is_perfect_square(y2))  { //y2 is a perfect square
	
	  int z2 = m -x2 - y2;
	
	  if (Is_perfect_square(z2) && z2<=y2) { // z2 is a perfect square
	      if(z2==0 && y2==0) deg_val += 6;
	      else if(z2==0 && y2!= 0) {
		if(x2==y2) deg_val += 3*4;
		else deg_val+=3*4*2;
	      
	      }
	      else {  // x, y, z > 0
		if(z2==y2 && y2==x2) deg_val += 2*2*2;
		else if(z2==y2 && y2 != x2) deg_val += 3*2*4;
		else if(x2==y2 && z2 != y2) deg_val += 3*2*4;
		else deg_val += 6*2*2*2;
	      
	      }
	  }
	}
      }
    }
  }
  
  
  
  return deg_val;
  
}


//Integrator: rectangular, trapezoidal, 1/3 Simpson, 3/8 Simpson
double w(int t, int Simps_ord) {

  if(Simps_ord == 1) return 1.0;

  else if(Simps_ord == 2) {

    if(t == 0) return 0.5;
    else return 1.0;

  }

  else if(Simps_ord == 3) {

    if(t==0) return 2.0/6.0;
    if(t%2==0) return 4.0/6.0;
    else return 8.0/6.0;

  }

  else if(Simps_ord == 4) {
    if(t==0) return 3.0/8.0;
    if (t%3 == 0) return 6.0/8.0;
    else return 9.0/8.0;
  }
  
  else crash("Simpson integrator with order: "+to_string(Simps_ord)+" not found");

  exit(-1);
  return 0.0;

}



/// Implements the trap to debug
void debug_loop()
{
  volatile int flag=0;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0) {
  printf("Entering debug loop on rank %d, flag has address %p please type:\n"
	 "$ gdb -p %d\n"
	 "$ set flag=1\n"
	 "$ continue\n",
	 rank,
	 &flag,
	 getpid());
  }
  
  if(rank==0)
    while(flag==0);
  
 MPI_Barrier(MPI_COMM_WORLD);

}


double Get_2l_alpha_s( double Q, double Nf, double Lambda) { //2loops

  double B0= 11.0 -(2.0/3.0)*Nf;
  double B1= 102.0 - (38.0/3)*Nf;
  double Q2= pow(Q,2);
  double L2= pow(Lambda,2);
  return (4*M_PI/(B0*log(Q2/L2)))*( 1.0 - (B1/pow(B0,2))*log( fabs(log(Q2/L2)) )/log(Q2/L2));
  
}


double Get_3l_alpha_s( double Q, double Nf, double Lambda) { //3loop alpha 

 
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= 0.0;
  double B4= 0.0;
  double t=log(Q*Q/(Lambda*Lambda));
  double l=log(t);
 
  double ret = (1.0/(B0*t))*( 1 - B1*l/(B0*B0*t) + (pow(B1,2)*(pow(l,2) -l -1.0) + B0*B2)/(pow(B0,4)*pow(t,2)));
  
  ret = ret + (1.0/(B0*t))*(  (pow(B1,3)*(-2*pow(l,3)+5*pow(l,2)+ 4.0*l -1) - 6.0*B0*B2*B1*l + B0*B0*B3)/(2.0*pow(B0,6)*pow(t,3)));
  
  ret = ret + (1.0/(B0*t))*( (18.0*B0*B2*pow(B1,2)*(2*l*l -l -1.0) + pow(B1,4)*(6*pow(l,4) -26*pow(l,3) -9*l*l + 24*l + 7.0))/(6.0*pow(B0,8)*pow(t,4)));
  
  ret = ret + (1.0/(B0*t))*( ( -pow(B0,2)*B3*B1*(12*l+1) + 2*pow(B0,2)*5*pow(B2,2))/(6*pow(B0,8)*pow(t,4)));

  ret = ret + (1.0/(B0*t))*( 2*pow(B0,2)*B0*B4)/(6*pow(B0,8)*pow(t,4));

  return ret;

}


double Get_4l_alpha_s( double Q, double Nf, double Lambda) { //5loop alpha 

  double zeta_3 = riemann_zeta(3);
  double zeta_4 = riemann_zeta(4);
  double zeta_5 = riemann_zeta(5);
 

  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double B4= (1.0/pow(4*M_PI,5))*( 8157455.0/16 + 621885.0*zeta_3/2.0 - 88209.0*zeta_4/2 -288090*zeta_5 + Nf*( -336460813.0/1944 - 4811164*zeta_3/81 + 33935*zeta_4/6 + 1358995*zeta_5/27) + pow(Nf,2)*( 25960913.0/1944 + 698531.0*zeta_3/81 - 10526*zeta_4/9 - 381760.0*zeta_5/81) + pow(Nf,3)*( -630559.9/5832 - 48722.0*zeta_3/243 + 1618.0*zeta_4/27 + 460.0*zeta_5/9) + pow(Nf,4)*( 1205.0/2916 -152.0*zeta_3/81))   ;
  double t=log(Q*Q/(Lambda*Lambda));
  double l=log(t);
 
  double ret = (1.0/(B0*t))*( 1 - B1*l/(B0*B0*t) + (pow(B1,2)*(pow(l,2) -l -1.0) + B0*B2)/(pow(B0,4)*pow(t,2)));
  
  ret = ret + (1.0/(B0*t))*(  (pow(B1,3)*(-2*pow(l,3)+5*pow(l,2)+ 4.0*l -1) - 6.0*B0*B2*B1*l + B0*B0*B3)/(2.0*pow(B0,6)*pow(t,3)));
  
  ret = ret + (1.0/(B0*t))*( (18.0*B0*B2*pow(B1,2)*(2*l*l -l -1.0) + pow(B1,4)*(6*pow(l,4) -26*pow(l,3) -9*l*l + 24*l + 7.0))/(6.0*pow(B0,8)*pow(t,4)));
  
  ret = ret + (1.0/(B0*t))*( ( -pow(B0,2)*B3*B1*(12*l+1) + 2*pow(B0,2)*5*pow(B2,2))/(6*pow(B0,8)*pow(t,4)));

  ret = ret + (1.0/(B0*t))*( 2*pow(B0,2)*B0*B4)/(6*pow(B0,8)*pow(t,4));

  return ret;

}

void Print_4l_alpha_s() {

  double Qmax= 100;
  double Qmin=0.1;
  double Qstep=0.01;
  double Q=Qmin;
  ofstream Print_alpha("../data/gm2/alphas_Q.dat");
  Print_alpha<<"#Q[GeV]  Nf=3    Nf=4"<<endl;
  while( Q <= Qmax) {
    Print_alpha<<Q <<" "<<Get_4l_alpha_s(Q,3)<<" "<<Get_4l_alpha_s(Q, 4)<<endl;
    Q+=Qstep;
  }
  Print_alpha.close();
  Q=Qmin;
  ofstream Print_alpha_t("../data/gm2/alphas_t_s_pi.dat");
  Print_alpha_t<<"#t[fm]  Nf=3    Nf=4"<<endl;
  while( Q <= Qmax) {
    double flav3= Get_4l_alpha_s(Q,3);
    double flav4= Get_4l_alpha_s(Q,4);
    Print_alpha_t<<(M_PI/Q)*0.197327 <<" "<<flav3<<" "<<(1.0 - exp(-pow( flav3/M_PI,2)))<<" "<<flav4<<" "<<(1.0 - exp( - pow(flav4/M_PI,2)))<<endl;
    Q+=Qstep;
  }

  Print_alpha_t.close();
  Q=Qmin;
  ofstream Print_alpha_t_s1("../data/gm2/alphas_t_s_1.dat");
  Print_alpha_t_s1<<"#t[fm]  Nf=3    Nf=4"<<endl;
  while( Q <= Qmax) {
    double flav3= Get_4l_alpha_s(Q,3);
    double flav4= Get_4l_alpha_s(Q,4);
    Print_alpha_t_s1<<(1.0/Q)*0.197327 <<" "<<flav3<<" "<<(1.0 - exp(-pow( flav3/1.0,2)))<<" "<<flav4<<" "<<(1.0 - exp( - pow(flav4/1.0,2)))<<endl;
    Q+=Qstep;
  }

  Print_alpha_t_s1.close();
  return;
}


double MS_bar_to_pole_mass( double mu, double Nf, double Lambda, double m, double mc ) {

  double d= 12.0/(33.0 - 2*Nf);
  double b1= -107.0*d*d/16.0 + 19.0*d/4.0;
  double b2 = -37117.0*pow(d,3)/1536.0 + 243.0*pow(d,2)/32.0 + 325.0*pow(d,2)/96.0;
  double c1 = 23.0*pow(d,2)/12.0 + (5.0/6.0 -b1)*d;
  double zeta_3 = riemann_zeta(3);
  double c2= (55.0*zeta_3/8.0 -2591.0/576.0)*pow(d,3) - (5.0*zeta_3/2.0 - 313.0/48.0)*pow(d,2) - 0.5*( 35.0/36.0 + b2)*d + 0.5*c1*(c1-b1);
  
  double K1= 4.0/3.0;
  double K2= 16.11 -1.04*(Nf-1) -1.04*(1.0 -m_MS_bar_m(mu,Nf, Lambda, mc)/m_MS_bar_m(mu,Nf,Lambda,m)); 
  //double K2= (M_PI*M_PI*log(2)/9.0) + 7.0*M_PI*M_PI/18.0 - zeta_3/6.0 + 3673.0/288.0 - (M_PI*M_PI/18.0 +71.0/144.0)*(Nf-1);

  //running of MS-bar quark mass
  auto C = [&c1,&c2,&d, &Lambda, &Nf](double x) {
    double y = Get_4l_alpha_s(x,Nf,Lambda)/(M_PI*d);
    return pow(y,d) + c1*pow(y,1+d) + c2*pow(y,2+d);
  };

  
  
  //ratio between pole mass M and MS-bar quark mass at scale M
  auto D = [&K1, &K2, &Lambda, &Nf](double x) {
    double y = Get_4l_alpha_s(x,Nf,Lambda)/(M_PI);
    return 1.0 + K1*y + K2*pow(y,2);
  };

  //solve equation r= M/(C(M)*D(M));
  //double r= m/C(mu);
 
  auto F= [&](double x) { return x - MS_bar_mass_evolutor(x,mu,Nf,Lambda,1)*D(x)*m ;};
  //cout<<"r : "<<r<<" F(xmin): "<<F(0.7*m)<<" F(xmax): "<<F(10*m)<<endl;

  return R_brent(F, 0.7*m, 5*m);
}

double pole_mass_to_MS_bar(double mu, double Nf, double Lambda, double Mpole) {

  double d= 12.0/(33.0 - 2*Nf);
  double b1= -107.0*d*d/16.0 + 19.0*d/4.0;
  double b2 = -37117.0*pow(d,3)/1536.0 + 243.0*pow(d,2)/32.0 + 325.0*pow(d,2)/96.0;
  double c1 = 23.0*pow(d,2)/12.0 + (5.0/6.0 -b1)*d;
  double zeta_3 = riemann_zeta(3);
  double c2= (55.0*zeta_3/8.0 -2591.0/576.0)*pow(d,3) - (5.0*zeta_3/2.0 - 313.0/48.0)*pow(d,2) - 0.5*( 35.0/36.0 + b2)*d + 0.5*c1*(c1-b1);
  
  double K1= 4.0/3.0;
  double K2= (M_PI*M_PI*log(2)/9.0) + 7.0*M_PI*M_PI/18.0 - zeta_3/6.0 + 3673.0/288.0 - (M_PI*M_PI/18.0 +71.0/144.0)*(Nf-1);

  //running of MS-bar quark mass
  auto C = [&c1,&c2,&d, &Lambda, &Nf](double x) {
    double y = Get_2l_alpha_s(x,Nf,Lambda)/(M_PI*d);
    return pow(y,d) + c1*pow(y,1+d) + c2*pow(y,2+d);
  };
  
  //ratio between pole mass M and MS-bar quark mass at scale M
  auto D = [&K1, &K2, &Lambda, &Nf](double x) {
    double y = Get_2l_alpha_s(x,Nf,Lambda)/(M_PI);
    return 1.0 + K1*y + K2*pow(y,2);
  };
  
  return Mpole/(D(Mpole)*MS_bar_mass_evolutor(Mpole,mu,Nf,Lambda,1));
  
  // double r= Mpole/(C(Mpole)*D(Mpole));

  //return C(mu)*r;
}


double m_MS_bar_m( double mu, double Nf, double Lambda, double m ) {

  //double c0 = 1.0/( (33.0-2.0*Nf)/(12.0) );  // pow(mu,c0) is LO evolutor  x^(1-c0) = m/mu^c0

  //double guess = pow( m/pow(mu,c0), 1.0/(1-c0));

  auto F = [&](double x) { return x - MS_bar_mass_evolutor(x,mu,Nf,Lambda,1)*m ;};

  double mm= R_brent(F, MS_bar_mass_evolutor(0.95*min(mu,m),mu, Nf,Lambda,1)*m, MS_bar_mass_evolutor(1.05*max(mu,m),mu, Nf,Lambda,1)*m );


  return  mm;
}


double MS_bar_to_pole_mass_bis( double mu, double Nf, double Lambda, double m , double mc) {

  double mm= m_MS_bar_m( mu,Nf, Lambda,m);
  double mcmc= m_MS_bar_m( mu,Nf, Lambda,mc);


  double A= Get_4l_alpha_s(mm, Nf, Lambda)/M_PI;

  //cout<<"mm: "<<mm<<" mcmc: "<<mcmc<<" A: "<<A<<endl;
  return mm*( 1 + A*(4.0/3) + pow( A, 2)*( 13.4434 - 1.0414*(Nf-1) -1.0414*(1- 4.0*(mcmc/mm)/3.0 )) + pow(A,3)*( 0.6527*pow(Nf,2) -26.655*Nf + 190.595));  

  
}


double pole_mass_to_MS_bar_bis(double mu, double Nf, double Lambda, double Mpole, double mc) {

  double mcmc= m_MS_bar_m( mu,Nf, Lambda,mc);

  auto F = [&](double x) {
   double  A= Get_4l_alpha_s(x, Nf, Lambda)/M_PI;
   return (x/Mpole)*(1 + A*(4.0/3) + pow( A, 2)*( 13.4434 - 1.0414*(Nf-1) -1.0414*(1- 4.0*(mcmc/x)/3.0 )) + pow(A,3)*( 0.6527*pow(Nf,2) -26.655*Nf + 190.595))  -1 ;
  };

  double mm = R_brent(F, 0.7*Mpole, 1.5*Mpole);
  
  return mm*MS_bar_mass_evolutor(mu,mm,Nf,Lambda,1);  //C(mu)/C(mm);
}


double MS_bar_to_MRS_mass(double mu, double Nf, double Lambda, double m, double mc) {


  double zeta_3 = riemann_zeta(3);
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b2= (B2/B0) -pow(B1/B0,2);
  double b3= 0.5*( (B3/B0) -pow(B1/B0,3));
  
  double R0=0.535;
  double mm= m_MS_bar_m( mu,Nf, Lambda,m);
  double mcmc= m_MS_bar_m( mu,Nf, Lambda,mc);
  double A= Get_4l_alpha_s(mm, 3, Get_Lambda_MS_bar(3));
  double Ag=A/( 1 + b2*pow(A,2) + b3*pow(A,3));

  Vfloat dr({-0.1106, -0.0340, 0.0966, 0.0162});
  double K=1;
  for (int n=0;n<(signed)dr.size();n++) K += pow(Ag,n+1)*dr[n];

  double z= mcmc/mm;

  double Mc_A1 = mcmc*(1.596 - 0.6285*z + 0.1777*pow(z,2));
  double Mc_A2 = -mm*( 1.0414 + 0.4444*log(z));
  double Mc_B1 = mcmc*(19.987 + 2.824*z + 1.288*pow(z,2) - (13.644+2.788*z-0.0343*pow(z,2))*log(z));
  double Mc_B2 = -mm*( 22.312 + (8.227 + 1.064*z -0.419*pow(z,2) +0.118*pow(z,3) -0.148*log(z))*log(z));

  double DMc= (Mc_A1+Mc_A2)*pow(A/M_PI,2) + (Mc_B1+Mc_B2)*pow(A/M_PI,3);

  
  return (K+R0*J_MRS(3,mm))*mm + DMc ;
}


pair<double, double> MS_bar_to_mm_and_MRS_mass(double mu, double Nf, double Lambda, double m, double mcmc) {


  double zeta_3 = riemann_zeta(3);
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b2= (B2/B0) -pow(B1/B0,2);
  double b3= 0.5*( (B3/B0) -pow(B1/B0,3));
  
  double R0=0.535;
  double mm= m_MS_bar_m( mu,Nf, Lambda,m);
  
  double A= Get_4l_alpha_s(mm, 3, Get_Lambda_MS_bar(3));
  double Ag=A/( 1 + b2*pow(A,2) + b3*pow(A,3));

  Vfloat dr({-0.1106, -0.0340, 0.0966, 0.0162});
  double K=1;
  for (int n=0;n<(signed)dr.size();n++) K += pow(Ag,n+1)*dr[n];

  double z= mcmc/mm;

  double Mc_A1 = mcmc*(1.596 - 0.6285*z + 0.1777*pow(z,2));
  double Mc_A2 = -mm*( 1.0414 + 0.4444*log(z));
  double Mc_B1 = mcmc*(19.987 + 2.824*z + 1.288*pow(z,2) - (13.644+2.788*z-0.0343*pow(z,2))*log(z));
  double Mc_B2 = -mm*( 22.312 + (8.227 + 1.064*z -0.419*pow(z,2) +0.118*pow(z,3) -0.148*log(z))*log(z));

  double DMc= (Mc_A1+Mc_A2)*pow(A/M_PI,2) + (Mc_B1+Mc_B2)*pow(A/M_PI,3);

   
  return make_pair( mm, (K+R0*J_MRS(3,mm))*mm + DMc) ;
}






double MRS_mass_to_MS_bar_mass(double mu, double Nf, double Lambda, double M_MRS, double mc) {

  double R0=0.535;
  double mcmc= m_MS_bar_m( mu,Nf, Lambda,mc);

  double zeta_3 = riemann_zeta(3);
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b2= (B2/B0) -pow(B1/B0,2);
  double b3= 0.5*( (B3/B0) -pow(B1/B0,3));
  

  auto F = [&](double mm) {
  double A= Get_4l_alpha_s(mm, 3, Get_Lambda_MS_bar(3));
  double Ag=A/( 1 + b2*pow(A,2) + b3*pow(A,3));

  Vfloat dr({-0.1106, -0.0340, 0.0966, 0.0162});
  double K=1;
  for (int n=0;n<(signed)dr.size();n++) K += pow(Ag,n+1)*dr[n];

  double z= mcmc/mm;

  double Mc_A1 = mcmc*(1.596 - 0.6285*z + 0.1777*pow(z,2));
  double Mc_A2 = -mm*( 1.0414 + 0.4444*log(z));
  double Mc_B1 = mcmc*(19.987 + 2.824*z + 1.288*pow(z,2) - (13.644+2.788*z-0.0343*pow(z,2))*log(z));
  double Mc_B2 = -mm*( 22.312 + (8.227 + 1.064*z -0.419*pow(z,2) +0.118*pow(z,3) -0.148*log(z))*log(z));

  double DMc= (Mc_A1+Mc_A2)*pow(A/M_PI,2) + (Mc_B1+Mc_B2)*pow(A/M_PI,3);


  return (K+R0*J_MRS(3,mm))*mm + DMc - M_MRS ;
  };

 double mm = R_brent(F,0.7*M_MRS, 1.3*M_MRS);
 return mm*MS_bar_mass_evolutor( mu, mm, 4, Lambda, 1);
}

double MRS_mass_to_mm(double mu, double Nf, double Lambda, double M_MRS, double mc) {

  double R0=0.535;
  double mcmc= m_MS_bar_m( mu,Nf, Lambda,mc);

  double zeta_3 = riemann_zeta(3);
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b2= (B2/B0) -pow(B1/B0,2);
  double b3= 0.5*( (B3/B0) -pow(B1/B0,3));
  

  auto F = [&](double mm) {
  double A= Get_4l_alpha_s(mm, 3, Get_Lambda_MS_bar(3));
  double Ag=A/( 1 + b2*pow(A,2) + b3*pow(A,3));

  Vfloat dr({-0.1106, -0.0340, 0.0966, 0.0162});
  double K=1;
  for (int n=0;n<(signed)dr.size();n++) K += pow(Ag,n+1)*dr[n];

  double z= mcmc/mm;

  double Mc_A1 = mcmc*(1.596 - 0.6285*z + 0.1777*pow(z,2));
  double Mc_A2 = -mm*( 1.0414 + 0.4444*log(z));
  double Mc_B1 = mcmc*(19.987 + 2.824*z + 1.288*pow(z,2) - (13.644+2.788*z-0.0343*pow(z,2))*log(z));
  double Mc_B2 = -mm*( 22.312 + (8.227 + 1.064*z -0.419*pow(z,2) +0.118*pow(z,3) -0.148*log(z))*log(z));

  double DMc= (Mc_A1+Mc_A2)*pow(A/M_PI,2) + (Mc_B1+Mc_B2)*pow(A/M_PI,3);


  return (K+R0*J_MRS(3,mm))*mm + DMc - M_MRS ;
  };

 double mm = R_brent(F,0.7*M_MRS, 1.3*M_MRS);
 return mm;

}

double J_MRS(double Nf, double mm) {

  double zeta_3 = riemann_zeta(3);
  
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*zeta_3) -(1078361.0/162.0 + 6508.0*zeta_3/27.0)*Nf + (50065.0/162.0 + 6472.0*zeta_3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));
  double b2= (B2/B0) -pow(B1/B0,2);
  double b3= 0.5*( (B3/B0) -pow(B1/B0,3));
  double A= Get_4l_alpha_s(mm, 3, Get_Lambda_MS_bar(3));
  double Ag= A/( 1 + b2*pow(A,2) + b3*pow(A,3));
  double y=  1.0/(2*B0*Ag);
  double M= exp( -y )/(2*B0);
  double b= B1/(2*B0*B0);

  auto Fn = [&y, &b](double n) { return (y/n)*(n-1-b)/(n-b);};

  bool converged=false;
  double sum= (-1/b);
  double F= sum;
  double prec_sum= sum;
  int n=1;
  while(!converged) {
    F *= Fn(n);
    sum += F;
    converged= (sum==prec_sum);
    prec_sum= sum;
    n++;
  }

  return sum*M;

}



double MS_bar_mass_evolutor( double mu1, double mu2, double Nf, double Lambda, int mode) { // computes m(mu1)/mu(mu2) 

  double z3= riemann_zeta(3);
  double z4= riemann_zeta(4);
  double z5= riemann_zeta(5);

  double g1= 4;
  double g2 = 202.0/3.0 -20.0*Nf/9.0;
  double g3 = 1249.0 + (-2216.0/27.0 - 160.0*z3/3.0)*Nf -140.0*pow(Nf,2)/81.0;
  double g4= 4603055.0/162.0 + 135680.0*z3/27.0 -8800.0*z5;
  g4= g4+ (-91723.0/27.0 - 34192.0*z3/9.0 + 880.0*z4 + 18400.0*z5/9.0)*Nf;
  g4= g4+ (5242.0/243.0 + 800.0*z3/9.0 - 160.0*z4/3.0)*Nf*Nf;
  g4= g4+ (-332.0/243.0 + 64.0*z3/27.0)*pow(Nf,3);

  g1 /=  4;
  g2 /=pow(4,2);
  g3 /= pow(4,3);
  g4 /= pow(4,4);

  //definition mu^2 * dm(mu)/dmu^2 = - g*m(mu)  -> dlog(m)/dlog(mu) = 2g,     g= Gi*[ alpha(mu)/pi ]^i

  //running of alpha  dalpha/dlog(mu) = -2B,   B= alpha*( Bi*alpha^i )   -> d(alpha/pi)/dlogmu = -2B' , B' = (alpha/p') * ( Bi*(pi)^i )*(alpha/pi)^i
  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
  double B1= (153.0 - 19*Nf)/(24.0*M_PI*M_PI);
  double B2= (2857.0 -5033.0*Nf/9.0 + 325.0*Nf*Nf/27.0)/(128.0*pow(M_PI,3));
  double B3= ( (149753.0/6.0 + 3564*z3) -(1078361.0/162.0 + 6508.0*z3/27.0)*Nf + (50065.0/162.0 + 6472.0*z3/81.0)*Nf*Nf + 1093.0*pow(Nf,3)/729.0)/(256.0*pow(M_PI,4));

  B0 *= M_PI;
  B1 *= pow(M_PI,2);
  B2 *= pow(M_PI,3);
  B3 *= pow(M_PI,4);
  


 

  //integrate the evolutor

  double a= Get_4l_alpha_s(mu2,Nf, Lambda)/M_PI;
  double b= Get_4l_alpha_s(mu1, Nf, Lambda)/M_PI;


  //return approximate solution if mode==1 : https://arxiv.org/pdf/hep-ph/9708255.pdf [Eq. 59]
  double c0=g1/B0;
  double c1=g2/B0;
  double c2=g3/B0;
  double c3=g4/B0;

  double b1=B1/B0;
  double b2=B2/B0;
  double b3=B3/B0;


  auto C_approx = [&](double x) {

    return pow(x,c0)*( 1 + (c1-b1*c0)*x + 0.5*pow(x,2)*( pow(c1-b1*c0,2) + c2 - b1*c1 + b1*b1*c0 -b2*c0) 
		       + ( (1.0/6)*pow(c1-b1*c0,3)+ 0.5*(c1-b1*c0)*(c2-b1*c1 + b1*b1*c0 -b2*c0) 
		       + (1.0/3)*(c3-b1*c2 +b1*b1*c1 - b2*c1 -pow(b1,3)*c0 + 2*b1*b2*c0 - b3*c0))*pow(x,3));
  };


  if(mode==1) return C_approx(b)/C_approx(a);

  //evolutor is \int_a^b  g/B d(alpha/pi)    a= alphas(mu2)/pi, b=alphas(mu1)/pi

  auto F= [&](double x) { return (g1*x + g2*pow(x,2)+ g3*pow(x,3) + g4*pow(x,4))/(x*( B0*x + B1*pow(x,2) + B2*pow(x,3) + B3*pow(x,4)));};  // F = g/B
    
  gsl_function_pp<decltype(F)> Fgsl(F);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function *G = static_cast<gsl_function*>(&Fgsl);
  double val, err;
  gsl_integration_qags(G, a,b , 0.0, 1e-4, 10000, w, &val, &err);
  gsl_integration_workspace_free (w);
  
  if(err/val > 1e-3) crash("in MS_bar_mass_evolutor, target precision 1e-6 not reached while integrating RG equations");

  return exp(val);
}


double evolutor_ZT_MS_bar( double mu1, double mu2, int nloops)  { //evolves the Z_T RC from Z_T(mu1) to Z_T(mu2) in the MSbar scheme at three loops, namely the function returns Z_T(mu2)/Z_T(mu1)

  


  double Lambda_QCD_MS_Bar_4= Get_Lambda_MS_bar(4);
  double Lambda_QCD_MS_Bar_5= Get_Lambda_MS_bar(5);
  double Lambda_QCD_MS_Bar_3= Get_Lambda_MS_bar(3);
  
  double CF= 4.0/3;
  double CA= 3;
  double TF= 0.5;
  double zeta_3=  riemann_zeta(3);
  double zeta_4=  riemann_zeta(4);
  double zeta_5=  riemann_zeta(5);

  
  
  auto Int= [&CF, &CA, &TF, &zeta_3, &zeta_4, &zeta_5, &Lambda_QCD_MS_Bar_3,  &Lambda_QCD_MS_Bar_4, &Lambda_QCD_MS_Bar_5, &nloops ](double mu) {

    auto alphas = [&Lambda_QCD_MS_Bar_3, &Lambda_QCD_MS_Bar_4, &Lambda_QCD_MS_Bar_5](double x) {
      if(x < 1.28) { return Get_4l_alpha_s(x, 3, Lambda_QCD_MS_Bar_3);}
      if(x < 4.18)  {return  Get_4l_alpha_s(x, 4, Lambda_QCD_MS_Bar_4);}

      return Get_4l_alpha_s(x,5,Lambda_QCD_MS_Bar_5);
    };
    
    double a =  alphas(mu)/(4.0*M_PI);
       
    int Nf=4;
    //if( mu > 1.28) Nf=4;
    //if(mu > 4.18) Nf=5;

    double anomalous_dim=0;
    if(nloops==1) anomalous_dim= (4.0/3.0)*a;
    else if(nloops==2) anomalous_dim=  (4.0/3.0)*a  - 2.0*( 26*Nf - 543)*a*a/27.0;
    else if(nloops==3)  anomalous_dim = (4.0/3.0)*a  - 2.0*( 26*Nf - 543)*a*a/27.0 - pow(a,3)*(36*pow(Nf,2) + 1440*zeta_3*Nf + 5240*Nf + 2784*zeta_3 - 52555)/81.0;
    else if(nloops==4) anomalous_dim= (4.0/3.0)*a  - 2.0*( 26*Nf - 543)*a*a/27.0 - pow(a,3)*(36*pow(Nf,2) + 1440*zeta_3*Nf + 5240*Nf + 2784*zeta_3 - 52555)/81.0 + (pow(a,4)/1458.0)*( 1152*zeta_3*pow(Nf,3) + 168.0*pow(Nf,3) + 66240*zeta_3*pow(Nf,2) -25920*zeta_4*pow(Nf,2) + 39844*pow(Nf,2) -1821984*zeta_3*Nf + 377568*zeta_4*Nf + 993600*zeta_5*Nf -3074758*Nf -742368*zeta_3 + 826848*zeta_4 -4018560*zeta_5 + 19876653);
    else crash("nloops must be <= 3");
    
    
    return -2*anomalous_dim/mu;
    
  };


  double val, err;

 
  double prec=1e-6;
  gsl_function_pp<decltype(Int)> integrand(Int);
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function *G = static_cast<gsl_function*>(&integrand);
  gsl_integration_qags(G, mu1, mu2,  0.0, prec, 1000, w, &val, &err);
  gsl_integration_workspace_free(w);
  if(fabs(err/val) > 5*prec) crash("In Compute_differential_decay_rate, cannot reach target precision: "+to_string_with_precision(prec,5));


  return exp(val);

  
		      

}

double Get_Lambda_MS_bar(int Nf) {

  double MZ= 91.1876;
  double alpha_at_MZ=0.1179;
  Vfloat thr({1.280, 4.203});

  assert(Nf >= 3);

  if(Nf == 5) {
  
  auto F = [&](double x) { return Get_4l_alpha_s(MZ, Nf, x) - alpha_at_MZ ;};

  double B0= (33.0-2.0*Nf)/(12.0*M_PI);
    
  double guess= MZ*exp(-1/(alpha_at_MZ*B0*2.0));
  
  return R_brent(F, 0.5*guess, 4*guess);

  }

  else {
    double mth= thr[Nf-3];
    double L_Nfp1_flav= Get_Lambda_MS_bar(Nf+1);
    //cout<<"L(" << (Nf+1) <<") = "<<L_Nfp1_flav<<endl;
    auto G = [&](double x) { return Get_4l_alpha_s(mth, Nf,x)*(1 - (11.0/(72*M_PI*M_PI))*pow(Get_4l_alpha_s(mth,Nf,x),2) ) - Get_4l_alpha_s(mth,Nf+1, L_Nfp1_flav);};
    return R_brent(G, 0.5*L_Nfp1_flav, 1.5*L_Nfp1_flav);
    
  }

  return 0.0;
}
  



