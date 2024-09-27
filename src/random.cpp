#include "../include/random.h"



using namespace std;


double gauss(Pfloat A, GaussianMersenne& B) {

  return A.first + A.second*B();
}

double gauss(double mean, double err,  GaussianMersenne& B) {

  return mean + err*B();
}


Vfloat Covariate(Eigen::MatrixXd Cov, Eigen::VectorXd Vec,  GaussianMersenne& B) {

  if(Cov.rows() != Cov.cols()) crash("In Covariate, Covariance matrix  is not a square matrix");
  Eigen::MatrixXd Cov_sqrt = Cov.sqrt();
  if(Cov.rows() != Vec.rows()) crash("In Covariate, sizes of covariance matrix and vector of averages do not match");

  Eigen::VectorXd G(Cov.rows());
  
  for(unsigned int i=0; i<G.rows();i++) G(i) = B();

  Eigen::VectorXd Res(Cov.rows());
  Res= Cov_sqrt*G;
  Vfloat T;

  for(unsigned int i=0; i < Vec.rows(); i++) T.push_back( Vec(i) + Res(i));

  return T;
}


Vfloat Covariate(Eigen::MatrixXd Cov,  GaussianMersenne& B) {

  if(Cov.rows() != Cov.cols()) crash("Covariance matrix in Covariate is not a square matrix");
  Eigen::MatrixXd Cov_sqrt = Cov.sqrt();

  Eigen::VectorXd G(Cov.rows());
  
  for(unsigned int i=0; i<G.rows();i++) G(i) = B();

  Eigen::VectorXd Res(Cov.rows());
  Res= Cov_sqrt*G;
  Vfloat T;

  for(unsigned int i=0; i < Cov.rows(); i++) T.push_back(  Res(i));

  return T;
}

