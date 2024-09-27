#ifndef __numerics__
#define __numerics__

#include "matplotlibcpp.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <map>
#include <string>
#include <functional>
#include <stdarg.h>
#include <numeric>
#include <sstream>
#include <cassert>
#include <utmpx.h>
#include <fenv.h>
#include <sched.h>
#include <mpi.h>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <bits/stdc++.h> 
#include <boost/algorithm/string.hpp> 
#include <ctime>
#include <chrono>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnMinimize.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/naive_monte_carlo.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <complex>


using namespace std;

typedef boost::minstd_rand base_generator_type;
typedef vector<int> Vint;
typedef vector<vector<int>> VVint;
typedef pair<int,int> Pint;
typedef pair<double,double> Pfloat;
typedef vector<pair<double,double>> VPfloat;
typedef vector<vector<pair<double,double>>> VVPfloat;
typedef vector<pair<int,int>> VPint;
typedef vector<double> Vfloat;
typedef vector<vector<double>> VVfloat;
typedef vector<long double> Vdouble;
typedef vector<vector<pair<int,int>>> VVPint;
typedef vector<vector<vector<vector<double>>>> VVVVfloat;
typedef vector<vector<vector<double>>> VVVfloat;
typedef vector<vector<vector<pair<double,double>>>> VVVPfloat;
typedef vector<vector<vector<vector<pair<double,double>>>>> VVVVPfloat;



//wrapper for lambdas with capture to gsl_monte_function
template< typename F >  class gsl_monte_function_pp : public gsl_monte_function {
public:
  gsl_monte_function_pp(const F& func) : _func(func) {
    f = &gsl_monte_function_pp::invoke;
    params=this;
    dim= 5;
  }
private:
  const F& _func;
  static double invoke(double x[], size_t dim, void *params) {
    vector<double> xv;
    for(int i=0; i<(signed)dim;i++) xv.push_back(x[i]);
    return static_cast<gsl_monte_function_pp*>(params)->_func(xv);
  }
};


//wrapper for lambdas with capture to gsl_function
template< typename F >  class gsl_function_pp : public gsl_function {
public:
  gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
  }
private:
  const F& _func;
  static double invoke(double x,void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
  }
};

//wrapper for lambdas with capture to gsl_function_fdf
template< typename F, typename F2, typename fdF >  class gsl_function_fdf_pp : public gsl_function_fdf {
public:
  gsl_function_fdf_pp(const F& func, const F2& df_func, const fdF& fdf_func) : _func(func), _df_func(df_func), _fdf_func(fdf_func) {
    f = &gsl_function_fdf_pp::invoke;
    df = &gsl_function_fdf_pp::invoke_df;
    fdf= &gsl_function_fdf_pp::invoke_fdf;
    params=this;
  }
private:
  const F& _func;
  const F2& _df_func;
  const fdF& _fdf_func;
  static double invoke(double x,void *params) {
    return static_cast<gsl_function_fdf_pp*>(params)->_func(x);
  }
   static double invoke_df(double x,void *params) {
    return static_cast<gsl_function_fdf_pp*>(params)->_df_func(x);
  }
  static void invoke_fdf(double x, void *params, double* f, double *df) {
    return static_cast<gsl_function_fdf_pp*>(params)->_fdf_func(x,f,df);
  }
};


const auto fake_func= [](const function<double(double)> &A) { return 0.0;};
const auto fake_func_d = [](double x) { return 0.0;};
void D(int k);
long long int ipow(int a,int n);
double fpow(double a, int n);
long long int fact(int n);
long long int BinomialCoeff(int n,int k);
double RatioPol(Vfloat &x_num, Vfloat &x_den, Vfloat &Num, Vfloat &Den, Vint &NumPow, Vint &DenPow, string MODE_POL);
void crash(string Message);
double eps(int k);
double F(double x, double a, int alpha);
void derivative(Vfloat &RES, Vfloat &INPUT, string MODE);
double FTAN(double x1, int t, int NT);
double FTAN_SYMM(double x1, int t, int NT);
double Root_Brent(double R, int nt, int NT);
double Root_Brent_sinh(double R, int nt, int NT);
double DoConstantFit(Vfloat &data, Vfloat &err);
double lin_interpolator(double y1, double y2, double Dx1, double Dx2,double Dx);
double quad_interpolator(double y1, double y2, double y3, double Dx1, double Dx2, double Dx3, double Dx);
void Print_To_File(const vector<string>& row_id, const VVfloat &data, string Path, string MODE, string Header);
double Get_4l_alpha_s(double Q, double Nf, double Lambda = 0.340);
double Get_3l_alpha_s(double Q, double Nf, double Lambda=0.340);
double Get_2l_alpha_s( double Q, double Nf, double Lambda=0.340);
void Print_4l_alpha_s();
double m_MS_bar_m( double mu, double Nf, double Lambda, double m );
double MS_bar_to_pole_mass(double mu, double Nf, double Lambda, double m, double mc);
double MS_bar_to_pole_mass_bis(double mu, double Nf, double Lambda, double m, double mc);
double pole_mass_to_MS_bar(double mu, double Nf, double Lambda, double Mpole);
double pole_mass_to_MS_bar_bis(double mu, double Nf, double Lambda,double Mpole, double mc);
double MS_bar_to_MRS_mass(double mu, double Nf, double Lambda, double m, double mc);
double MRS_mass_to_MS_bar_mass(double mu, double Nf, double Lambda, double M_MRS, double mc);
double MRS_mass_to_mm(double mu, double Nf, double Lambda, double M_MRS, double mc);
double MS_bar_mass_evolutor(double mu1, double mu2, double Nf, double Lambda,  int mode);
pair<double, double> MS_bar_to_mm_and_MRS_mass(double mu, double Nf, double Lambda, double m, double mc);
double J_MRS(double Nf, double mm) ;
double evolutor_ZT_MS_bar( double mu1, double mu2, int nloops);
double Get_Lambda_MS_bar(int Nf);

double w(int t, int Simps_ord);







template <typename T>
string to_string_with_precision(T a_value, const int n)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

template <typename T>
void printV(const vector<T> &A, string B, bool mode) {
  int size = A.size();
  cout.precision(10);
  cout << B << endl;
  if(mode ==0 )for(int i = 0;i < size; i++) cout<< A[i] << " ";
  else for(int i=0;i<size;i++) cout<<i<<"  "<<A[i]<<endl;
  cout << endl;
  return;

}

template <typename T>

T Compute_scalar_product( const vector<T> & A,const vector<T> &B) {

T res =0;

if(A.size() != B.size()) crash("In compute scalar product, size of the two vectors do not coincide");

for (unsigned int i=0; i < A.size(); i++) res += A[i]*B[i];

return res;

}

template <typename T>
vector<T> Multiply_vectors( const vector<T> & A, const vector<T> & B) {

  vector<T> res;
  if(A.size() != B.size()) crash("In Multiply_vectors, the size of the two vectors do not coincide");

  for(unsigned int i=0; i < A.size(); i++) res.push_back(A[i]*B[i]);

  return res;

}


template <typename T>
vector<T> Sum_vectors (const vector<T> & A,const vector<T> &B) {

 vector<T> res;
  if(A.size() != B.size()) crash("In Sum_vectors, the size of the two vectors do not coincide");

  for(unsigned int i=0; i < A.size(); i++) res.push_back(A[i]+B[i]);

  return res;


}

template <typename T>

vector<vector<T>> Sum_Vvectors( const vector<vector<T>> & A,const  vector<vector<T>> &B) {

  vector<vector<T>> res;

  if(A.size() != B.size()) crash("In Multiply_Vvectors, the size of the two vectors do not coincide");


  for(unsigned int i=0; i<A.size();i++) res.push_back( Sum_vectors(A[i], B[i]));

  return res;



}



template <typename T,
	  typename...V>
T summ_master(const T& t,const V&...v)
{
  return (v+...+t);
}

template <typename T,
	  typename...V>
auto summ_master(const std::vector<T>& t0,const std::vector<V>&...t)
{
   
  static_assert((std::is_same_v<T,V> and...),"Needs to call with the same container type");
  
    
  std::vector<T> res(t0.size());
  
  for(std::size_t i=0;i<t0.size();i++)
    res[i]=summ_master(t0[i],t[i]...);
  
  return res;
}


template <typename T,
	  typename...V>
T multiply_master(const T& t,const V&...v)
{
  return (v*...*t);
}

template <typename T,
	  typename...V>
auto multiply_master(const std::vector<T>& t0,const std::vector<V>&...t)
{
   
  static_assert((std::is_same_v<T,V> and...),"Needs to call with the same container type");
  
    
  std::vector<T> res(t0.size());
  
  for(std::size_t i=0;i<t0.size();i++)
    res[i]=multiply_master(t0[i],t[i]...);
  
  return res;
}




template <typename T>
double R_brent( T&& F, double xmin, double xmax) {

  double Precision = 1e-6;
  double delta = 0.01;
  //solve the equation F(X) = 0 using the Brent method!
  //initialize iteration
  double b=xmax;
  double a=xmin;
 
  if (F(a)*F(b) > 0) {crash("Initial conditions in R_brent not valid: (xmin,xmax) = ("+to_string_with_precision(xmin,5)+","+to_string_with_precision(xmax,5)+") , (f(xmin), f(xmax)) = ("+to_string_with_precision(F(xmin),5)+","+to_string_with_precision(F(xmax),5)+")");
  }
 
 



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



template <typename T>
vector<T> Multiply_vector_by_scalar(const vector<T> & A, T B) {

  vector<T> res;

  for(auto & el: A) res.push_back(el*B);

  return res;


}

template <typename T>

vector<vector<T>> Multiply_Vvector_by_scalar( const vector<vector<T>> &A, T B) {

  vector<vector<T>> res;

  for( auto & el: A) res.push_back ( Multiply_vector_by_scalar( el, B));

  return res;
}

template <typename T>
void Transpose_VV(vector<vector<T>> &A) {

  if(A.size() == 0) crash("Transpose_VV called with empty matrix.");
  

  vector<vector<T>> B;
  B.resize(A[0].size());
  for(auto &B_i:B) B_i.resize(A.size());

  for(unsigned int i=0; i<A.size(); i++) {
    if(A[i].size() != B.size() ) crash("Matrix A in Transpose_VV has an invalid size");
    for(unsigned int j=0; j<A[i].size(); j++) B[j][i] = A[i][j];
  }



  A=B;
  return;



}


template <typename T>
vector<T> slicing(const vector<T>& arr, int X, int Y) {
  
    // Starting and Ending iterators 
    auto start = arr.begin() + X; 
    auto end = arr.begin() + Y + 1; 
  
    // To store the sliced vector 
    vector<T> result(Y - X + 1); 
  
    // Copy vector using copy function() 
    copy(start, end, result.begin()); 
  
    // Return the final sliced vector 
    return result; 
}


template <typename T>
vector<T> external_prod( const vector<T> &arr1, const vector<T> &arr2) {

  if (arr1.size() != 3 || arr2.size() != 3) crash("Invalid call to external_prod, vectors sizes != 3");

  vector<T> res(3);

  res[0] = arr1[1]*arr2[2] - arr1[2]*arr2[1];
  res[1] = arr1[2]*arr2[0] - arr1[0]*arr2[2];
  res[2] = arr1[0]*arr2[1] - arr1[1]*arr2[0];

  return res;


};

template <typename T>
T Kahan_sum(const vector<T> &input) {

  T sum= 0.;
  T c = 0.;

  for (auto & val: input) {
    T y = val- c;
    T t = sum + y;
    c =(t-sum) -y;
    sum = t;
  }
   
  return sum;

};



void cascade_resize( vector<vector<double>>& arr, const Vint& A ); 
void cascade_resize( vector<vector<vector<double>>>& arr,const Vint &A);
void cascade_resize( vector<vector<vector<vector<double>>>>& arr,const Vint &A);
void cascade_resize( vector<vector<vector<vector<vector<double>>>>>& arr,const Vint &A);
bool Is_perfect_square(int x);
int degeneracy(int m);



//define special functions
const auto g1_l = [](double x) -> double {

		      double n_max= 50;

		      double res=0.0;

		      //double res_asympt=0.0;

		      //for(int n=1;n<=n_max;n++) res_asympt += 4.0*sqrt(M_PI/2.0)*degeneracy(n)*exp(-sqrt(n)*x)/pow(sqrt(n)*x,1.5);

		      //return res_asympt;

		      for(int n=1; n<=n_max;n++) {
			
			res += (4.0*degeneracy(n)/(sqrt(n*1.0)*x))*boost::math::cyl_bessel_k(1, sqrt(n*1.0)*x);
			
		      }
		     
		      return res;
		    };


const auto g2_l = [](double x) -> double {

		      double n_max= 50;

		      double res=0.0;

		      //double res_asympt=0.0;

		      //for(int n=1;n<=n_max;n++) res_asympt += 4.0*sqrt(M_PI/2.0)*degeneracy(n)*exp(-sqrt(n)*x)/pow(sqrt(n)*x,1.5);

		      //return res_asympt;

		      for(int n=1; n<=n_max;n++) {
			
			res += (4.0*degeneracy(n)/(pow(sqrt(n*1.0)*x,2)))*boost::math::cyl_bessel_k(2, sqrt(n*1.0)*x);
			
		      }
		     
		      return res;
		    };



//CDH FORMULAE FOR Mpi and fpi
const auto Cf1 = [](double l1, double l2, double l3, double l4) { return -(7.0/9.0) + 2.0*l1 + (4.0/3.0)*l2 - 3.0*l4;};
const auto Cf2 = [](double l1, double l2, double l3, double l4) { return  112.0/9.0 - (8.0/3.0)*l1 - (32.0/3.0)*l2;};
const auto Cm1 = [](double l1, double l2, double l3, double l4) { return (-55.0/18.0) + 4.0*l1 +(8.0/3.0)*l2 - (5.0/2.0)*l3 -2.0*l4;};
const auto Cm2 = [](double l1, double l2, double l3, double l4) {return (112.0/9.0) - (8.0/3.0)*l1 -(32.0/3.0)*l2;};
const auto Sf1 = [](double s0, double s1, double s2, double s3) { return (4.0/3.0)*s0 - (13.0/6.0)*s1;};
const auto Sf2 = [](double s0, double s1, double s2, double s3) { return -(40.0/3.0)*s0 +4.0*s1 +(8.0/3.0)*s2 + (13.0/3.0)*s3;};
const auto Sm1 = [](double s0, double s1, double s2, double s3) { return s0*13.0/3.0;};
const auto Sm2 = [](double s0, double s1, double s2, double s3) { return -(40.0/3.0)*s0 - (32.0/3.0)*s1 - (26.0/3.0)*s2;};

const auto Cf1_log = []() { return 2.0 + 4.0/3.0 - 3.0;};
const auto Cf2_log = []() { return -8.0/3.0 -32.0/3.0;};
const auto Cm1_log = []() { return 4.0 + 8.0/3.0 - 5.0/2.0 - 2.0;};
const auto Cm2_log = []() { return -8.0/3.0 -32.0/3.0;}; 
   

//debug
void debug_loop();

#endif


