#ifndef __kernel__
#define __kernel__

#include "numerics.h"

double logplus(double x, double y, double rl);

double logminus(double x, double y, double rl);

double ptrate(double x, double y, double rl, double rl2, double mk, double fk);

double kernV(double x, double y, double rl, double rl2, double mk, double fk);

double kernA(double x, double y, double rl, double rl2, double mk, double fk);

double kern1(double x, double y, double rl, double rl2, double mk, double fk);

double kern2(double x, double y, double rl, double rl2, double mk, double fk);

double kernVV(double x, double y, double rl, double rl2, double mk);

double kernAA(double x, double y, double rl, double rl2, double mk);

double kern11(double x, double y, double rl, double rl2, double mk);
  
double kern22(double x, double y, double rl, double rl2, double mk);

double kernA1(double x, double y, double rl, double rl2, double mk);

double kern12(double x, double y, double rl, double rl2, double mk);

long double Compute_square_amplitude_extended( double h1, double h2, double fA, double fV, double he1, double he2, double feA, double feV, double MK, double fk,double rl,double xk, double xq, double A, double B, double y, string MODE);

long double Compute_square_amplitude_different_lepton(double h1, double h2, double fA, double fV, double MK, double fk,double rl, double rll, double xk, double xq, double a, double b, double Y, string MODE );

long double Compute_square_amplitude_extended_v2( double h1, double h2, double fA, double fV, double he1, double he2, double feA, double feV, double MK, double fk,double rl,double xk, double xq, double A, double B, double y, string MODE);

long double Compute_square_amplitude_different_lepton_v2(double h1, double h2, double fA, double fV, double MK, double fk,double rl, double rll, double xk, double xq, double a, double b, double Y, string MODE );


















#endif
