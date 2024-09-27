#ifndef __pt3_momenta__
#define __pt3_momenta__



#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "header_file_virph.h"
#include "LatInfo.h"
#include "T_min.h"
#include "input.h"

using namespace std;



class pt3_momenta {

 public:
  pt3_momenta() {}
 pt3_momenta(Vfloat A, Vfloat B, Vfloat C) : theta_0(A), theta_s(B), theta_t(C) {
    if (A.size() != 3 || B.size() != 3 || C.size() != 3) crash("In class pt3_momenta, call to constructor is invalid: cannot build 3_momenta using vectors of size != 3");
  }
  pt3_momenta(double A, double B, double C, double mu1, double mu2, double virt, double L , double T) {
    theta_0 = {0.0, 0.0, A};
    theta_s = {0.0, 0.0, B};
    theta_t = {0.0, 0.0, C};
    this->mu1 = mu1;
    this->mu2 = mu2;
    this->virt_val = virt;
    this->l = L;
    this->nt = T;
  }
 pt3_momenta(Vfloat A, Vfloat B, Vfloat C, double mu1, double mu2, double virt, double L,double T): theta_0(A), theta_s(B), theta_t(C) {
    if(A.size() != 3 || B.size() != 3 || C.size() != 3) crash("In class pt3_momenta, call to constructor is invalid: cannot build 3_momenta using vectors of size != 3");
    this->mu1= mu1;
    this->mu2 = mu2;
    this->virt_val = virt;
    this->l = L;
    this->nt = T;
  }


  Vfloat k() {
    Vfloat res({this->theta_0[0]-this->theta_t[0], this->theta_0[1]-this->theta_t[1], this->theta_0[2]-this->theta_t[2]});
    res= Multiply_vector_by_scalar(res, 2*M_PI/l);
    for(auto & k_i: res) k_i = 2.0*sin(k_i/2.0);
    return res;
  }
  Vfloat p() {
    Vfloat res({theta_0[0]-theta_s[0], theta_0[1]-theta_s[1], theta_0[2]-theta_s[2]});
    res= Multiply_vector_by_scalar(res, 2*M_PI/l);
    for(auto & p_i: res) p_i = 2.0*sin(p_i/2.0);
    return res;
  }
  double Egamma() {
    double e=0.0;
    Vfloat k=this->k();
    for(auto &k_x: k ) e += pow(k_x,2);
    double eg = 2.0*asinh(0.5*sqrt(e));
    return sqrt( eg*eg + pow(this->virt(),2));
    //return 2.0*asinh(sqrt(e + pow(2*sinh( this->virt()/2.0)  ,2) )/2.0);
  }

  double Egamma_asinh() {
     double e=0.0;
     Vfloat k=this->k();
     for(auto &k_x: k ) e += pow(k_x,2);
     return 2.0*asinh(sqrt(e + pow(2*sinh( this->virt()/2.0)  ,2) )/2.0);


  }


  double k_mod() {
    double mod=0;
    for(auto kx: k()) mod += pow(kx,2);
    return mod;
  }

  double p_mod() {
    double mod=0;
    for(auto px:p()) mod += pow(px,2);
    return mod;
  }
    
  distr_t x_gamma(distr_t& Mass) {
    distr_t xg(Mass.UseJack);
    double pk= Compute_scalar_product(this->k(),this->p());
    return 2.0*(this->Egamma()/Mass) -2.0*pk/(Mass*Mass);
  }

  distr_t x_gamma_off(distr_t& Mass) {
    distr_t xg(Mass.UseJack);
    double pk= Compute_scalar_product(this->k(),this->p());
    return 2.0*(this->Egamma()/Mass) -2.0*pk/(Mass*Mass) -this->virt()*this->virt()/(Mass*Mass);
  }

  distr_t E(distr_t& Mass) {
    distr_t e(Mass.UseJack);
    double p2=0;
    for(auto & p_x: p()) p2 += p_x*p_x;
    for(auto & el : Mass.distr) e.distr.push_back(sqrt( el*el + p2));
    return e;
  }

  string name() {
    return to_string_with_precision(theta_0[2],4)+"_"+to_string_with_precision(theta_s[2],4)+"_"+to_string_with_precision(theta_t[2],4)+"_"+to_string_with_precision( mu(0),4)+"_"+to_string_with_precision(mu(1),4)+"_"+to_string_with_precision(virt_val,4);
  }

  Vfloat Theta(int i) {
    if(i>=3) crash("In pt3_momenta::Theta(int i), index i is larger than 2. Exiting...");
    if(i==0) return theta_0;
    else if(i==1) return theta_s;
    else return theta_t;
  }

  double mu(int i) {
    if (i==0) return mu1;
    else if(i==1) return mu2;
    else crash("In pt3_momenta::mu(int i), index i is not 0 or 1. Exiting...");
    return 0;
  }

 

  double L() { return this->l;}

  double Nt() { return this->nt;}

  double virt() { return this->virt_val;}
    

  void Set_theta(Vfloat th, int i) {
    if(i==0) this->theta_0 = th;
    else if(i==1) this->theta_s = th;
    else if(i==2) this->theta_t = th;   
    else crash("In pt3_momenta::Set_theta(Vfloat, int i), index i != 0,1,2. Exiting...");
    return ;
  }

  void Set_theta(Vfloat th0, Vfloat ths, Vfloat tht) {
    if(th0.size() != 3 || ths.size() != 3 || tht.size() != 3) crash("In pt3_momenta::Set_theta(Vfloat, Vfloat, Vfloat) vectors are not 3_vectors");
    this->theta_0 = th0;
    this->theta_s = ths;
    this->theta_t = tht;
    return;
  }
  void Set_mu(double mu1, double mu2) {
    this->mu1 = mu1;
    this->mu2 = mu2;
    return;
  }


  void Set_virt(double k2) { this->virt_val = k2; return;}


  void Set_L(double length) { this->l = length; return;}

  void Set_Nt(double T) {this->nt = T; return;}
  
 private:
  Vfloat theta_0;
  Vfloat theta_s;
  Vfloat theta_t;
  double mu1, mu2, l;
  double virt_val;
  double nt;
};


class pt2_momenta {


 public:
  pt2_momenta() {}
  pt2_momenta(Vfloat th0, Vfloat ths, double m1, double m2, double L, double a,double b) : theta_0(th0), theta_s(ths), mu1(m1), mu2(m2), l(L), i0(a), is(b) {
    if(th0.size() != 3 || ths.size() != 3 ) crash("In class pt2_momenta, call to constructor is invalid. Trying to build three momenta with vectors of size != 3");
  }




  Vfloat Theta(int i) {
    if(i>=2) crash("In pt2_momenta::Theta(int i), index i is larger than 1. Exiting...");
    if(i==0) return theta_0;
    else return theta_s;
  }

  double mu(int i) {
    if (i==0) return mu1;
    else if(i==1) return mu2;
    else crash("In pt2_momenta::mu(int i), index i is not 0 or 1. Exiting...");
    return 0;
  }


  double L() { return this->l;}



  void Set_theta(Vfloat th0, Vfloat ths) {
    if(th0.size() != 3 || ths.size() != 3) crash("In pt2_momenta::Set_theta( Vfloat, Vfloat) vectors are not 3_vectors");
    this->theta_0 = th0;
    this->theta_s = ths;
    return;
  }
  void Set_mu(double mu1, double mu2) {
    this->mu1 = mu1;
    this->mu2 = mu2;
    return;
  }


  void Set_L(double L) { this->l = L; return;}

  bool Is_k0() {
    double sum=0;
    for(int i=0; i< 3; i++) sum += theta_0[i] + theta_s[i];
    if(sum< eps(15)) return true;
    else return false;
  }

  int Get_i0() { return i0;}
  int Get_is() { return is;}

  void Set_i0(int I0) { this->i0 = I0;}
  void Set_is(int IS) { this->is = IS;}
  
 private:
  Vfloat theta_0;
  Vfloat theta_s;
  double mu1, mu2, l, i0, is;
  

};


class pt2_momenta_list {
  
 public:
  pt2_momenta_list() {}
  pt2_momenta_list(int Nens) : mom(Nens) {}
  vector<vector<pt2_momenta>> mom;
  Vfloat Nt;
};


class pt3_momenta_list {

 public:
  pt3_momenta_list() {}
  pt3_momenta_list(int Nens) : mom(Nens) {}
  vector<vector<pt3_momenta>> mom;
  Vfloat Nt;
  


};


void Add_to_mom_list(pt3_momenta_list &M, struct header_virph &header, double& L); 

void Add_to_mom_list(pt2_momenta_list &M, struct header_virph &header, double& L);



#endif
