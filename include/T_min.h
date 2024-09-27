#ifndef __T_min__
#define __T_min__


#include "numerics.h"


using namespace std;


class fit_t_res {

 public:
  fit_t_res() {}

  Vfloat pars;
  double chi2;
  int Iterations;
  bool Status;
  

};


class T_fit : public ROOT::Minuit2::FCNBase {


public:
 T_fit(const Vfloat &X,const Vfloat &Y,const Vfloat &Y_err ):  X(X), Y(Y), Y_err(Y_err), theErrorDef(1.0), NumberOfMeasurements(X.size())  {}
 T_fit() : theErrorDef(1.0), NumberOfMeasurements(1) {}
 T_fit(const Vfloat &X,const Vfloat &Y):  X(X), Y(Y), theErrorDef(1.0), NumberOfMeasurements(X.size())  {for(int imeas=0;imeas<NumberOfMeasurements;imeas++) Y_err.push_back(1.0);}
 
 
  virtual ~T_fit() {}
  function<double(const Vfloat &ip, double x, int imeas)> ansatz;
  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const Vfloat& par) const;
  void add_pars(double a) {init_pars.push_back(a);}
  void add_pars(const Vfloat& a) { for (auto &el : a) init_pars.push_back(el);}
  void add_par_errs(double a) {init_par_errs.push_back(a);}
  void add_par_errs(const Vfloat& a) { for (auto &el : a) init_par_errs.push_back(el);}

  fit_t_res fit();


  Vfloat X;
  Vfloat Y;
  Vfloat Y_err;



private:

  
double theErrorDef;
double NumberOfMeasurements;
Vfloat init_pars;
Vfloat init_par_errs;





};


#endif

