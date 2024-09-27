#ifndef __ChPT_form_factors__
#define __ChPT_form_factors__



#include "numerics.h"
#include "random.h"


using namespace std;

double H(double xk2, double m, double mu);
void Compute_ChPT_form_factors( function<double(double,double)>& H1, function<double(double,double)>& H2, function<double(double,double)>& FA, function<double(double,double)>& FV);




#endif 
