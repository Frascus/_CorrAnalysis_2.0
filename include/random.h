#ifndef __random__
#define __random__

#include "numerics.h"



using namespace std;



class RandomMersenne {

 public:
 RandomMersenne(int seed, int size, int Length, int Nbranches) : gen(seed) ,rand(0,size) {
    this->Length=Length;
    gen_index.resize(Nbranches, 0);
    for(int i=0; i<this->Length*Nbranches;i++) RandFlow.push_back(rand(gen));
  }
 RandomMersenne() {}
 RandomMersenne(int seed, int size) : gen(seed), rand(0,size) {}
  
  int operator()(int ibranch) {
    if((signed)gen_index.size() < ibranch+1) crash("Random Mersenne called in stream mode without previous initialization");
    gen_index[ibranch]++;
    if((signed)RandFlow.size() == 0) crash("Random Mersenne called in stream mode without previous initialization");
    return RandFlow[Length*ibranch + gen_index[ibranch]-1];
  }
  int operator()() {
    boost::variate_generator<boost::mt19937&, boost::uniform_int<>> c(gen, rand);
    return c();
  }


 private:
  boost::mt19937 gen;
  boost::uniform_int<int> rand;
  Vint RandFlow;
  Vint gen_index;
  int Length;
};


class GaussianMersenne {

 public:
 GaussianMersenne(int seed, int Length, int Nbranches) : gen(seed) {
    this->Length=Length; 
    gen_index.resize(Nbranches, 0);
    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >    a(gen, boost::normal_distribution<>());
    for(int i=0; i<this->Length*Nbranches;i++) GaussFlow.push_back(a());
  }

 GaussianMersenne() {} 
 GaussianMersenne(int seed) : gen(seed)  {std_dev=1.0;}
 GaussianMersenne(int seed, double std_dev) : gen(seed), std_dev(std_dev) {}
    
  
  
  double operator()(int ibranch) {
    if((signed)gen_index.size() < ibranch+1) crash("Gaussian Mersenne called in stream mode without previous initialization");
    gen_index[ibranch]++;
    if( (signed)GaussFlow.size() == 0) crash("Gaussian Mersenne called in stream mode without previous initialization");
    return GaussFlow[Length*ibranch + gen_index[ibranch]-1];
  }
  double operator()() {
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<>> b(gen, boost::normal_distribution<>(0.0, std_dev));
    return b();
  }
  
  



 private:
  boost::mt19937 gen;
  Vfloat GaussFlow;
  Vint gen_index;
  int Length;
  double std_dev;
} ;


class UniformMersenne  {

public:

 UniformMersenne(int seed) : gen(seed), x(1), y(1)  {}
 UniformMersenne(int seed, double x, double y) : gen(seed), x(x), y(y) {}
    
  
  double operator()() {
    boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<>> b(gen, boost::random::uniform_real_distribution<>(x, y));
    return b();
  }
  
  
 private:
  boost::mt19937 gen;
  double x;
  double y;
  
};


double gauss(Pfloat A, GaussianMersenne& B);
double gauss(double mean, double err,  GaussianMersenne& B);


Vfloat Covariate(Eigen::MatrixXd Cov, Eigen::VectorXd Vec,  GaussianMersenne& B);
Vfloat Covariate(Eigen::MatrixXd Cov,  GaussianMersenne& B);























#endif
