#ifndef __stat__
#define __stat__

#include "numerics.h"
#include "random.h"


using namespace std;



Pfloat BootAve(const Vfloat &BootDistr);
double Compute_boot_cov(const Vfloat& A,const Vfloat& B);
Pfloat JackAve(const Vfloat &JackDistr);
double Compute_jack_cov(const Vfloat&,const Vfloat&);
void Compute_covariance_matrix(bool UseJack, Eigen::MatrixXd& Cov, int nargs,...);
void Compute_correlation_matrix(bool UseJack, Eigen::MatrixXd& Corr, int nargs,...);
void Compute_autocorrelation_time(const Vfloat& data, string Path, string Tag);
  





class distr_t {

 public:
  distr_t()  {UseJack = 1;}
  distr_t(bool a) : UseJack(a) {}
  distr_t(bool a, Vfloat b) : UseJack(a), distr(b) {}
  distr_t(bool a, int size) : UseJack(a), distr(size,0.0) {}

 //////////////////////////////////////////

 ////////////////////////////////////

 double ave() const;
 double err() const;
 Pfloat ave_err() const;
 int size() const;
 static distr_t f_of_distr(function<double(double x)> F,const distr_t& D) {
   distr_t res(D.UseJack);
   for(auto &d: D.distr) res.distr.push_back(F(d));

   return res;
 }
 

 bool UseJack;
 Vfloat distr;
};


//operator overloading

distr_t operator*(const distr_t& A, const distr_t& B);
 
distr_t operator+(const distr_t& A,const  distr_t& B);
 
distr_t operator-(const distr_t& A, const distr_t& B);

distr_t operator/(const distr_t& A, const distr_t& B);

double operator%(const distr_t& A,const  distr_t& B) ;

distr_t operator*(const distr_t& A, double B);

distr_t operator*(double B,  const distr_t& A);

distr_t operator+(const distr_t& A, double B);

distr_t operator+(double B,const distr_t& A);
distr_t operator-(const distr_t& A, double B) ;

distr_t operator-(double B,const distr_t& A);
distr_t operator/(const distr_t& A, double B) ;

distr_t operator/(double B,const distr_t& A);

distr_t operator*(const distr_t& A,const Vfloat& B);
 
distr_t operator+(const distr_t& A, const Vfloat& B) ;
 
distr_t operator-(const distr_t& A, const Vfloat& B) ;
distr_t operator/(const distr_t& A,const Vfloat& B) ;

distr_t operator*(const Vfloat &B,const distr_t& A) ;
distr_t operator+(const Vfloat &B, const distr_t& A);
distr_t operator-(const Vfloat& B, const distr_t& A) ;
distr_t operator/(const Vfloat& B, const distr_t& A);

// end operator overloading



class distr_t_list {

 public:
  distr_t_list() {UseJack=1;}
  distr_t_list( bool sampling_type, int size) :  UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) distr_list.emplace_back(sampling_type);
  }
  distr_t_list( bool sampling_type, VVfloat raw_distr_list) : UseJack(sampling_type) {
    for (unsigned int i_distr=0;i_distr<raw_distr_list.size(); i_distr++) this->distr_list.emplace_back(sampling_type, raw_distr_list[i_distr]);
  }
  distr_t_list(bool sampling_type) : UseJack(sampling_type) {}
  distr_t_list(int size, distr_t A) : UseJack(A.UseJack),  distr_list(size,A) {}
  distr_t_list(int sampling_type, int size, int sample_size) : UseJack(sampling_type) {
    for(int i_distr=0; i_distr<size;i_distr++) this->distr_list.emplace_back(sampling_type,sample_size);}
  

 distr_t operator[](int k) {
    if(k>= this->size()) crash("In distr_t_list::operator[](int k) k >= size");
    return this->distr_list[k];
  }

  //////////////////////////////////

 
 Vfloat ave() const ; 
 double ave(int i_distr) const;
 Vfloat err() const;
 double err(int i_distr) const;
 vector<Pfloat> ave_err() const;
 Pfloat ave_err(int i_distr) const ;
 int size() const;
 Vfloat Get_distr_index(int k) const;
 static distr_t_list f_of_distr_list(function<double(double a, double b)> F,const distr_t_list& D_List)  {
   distr_t_list RET(D_List.UseJack, D_List.size());
   for(int t=0; t < RET.size();t++) {
     for(int k=0;k< D_List.distr_list[t].size();k++) RET.distr_list[t].distr.push_back( F(D_List.distr_list[t].distr[k], t));
   }

   return RET;
 }
 static distr_t_list f_of_distr(function<double(double a, double b, double c)> F,const distr_t& D, int size)  {
   distr_t_list RET(D.UseJack, size);
   for(int t=0; t < RET.size();t++) {
     for(int k=0;k< D.size();k++) RET.distr_list[t].distr.push_back( F(D.distr[k], t, size));
   }

   return RET;
 }

  static distr_t_list derivative(const distr_t_list& D_List, int mode) {
    
    int NT= D_List.size();
    distr_t_list RET(D_List.UseJack, NT);
    for(int t=0; t<RET.size();t++) {
      for(int k=0;k< D_List.distr_list[t].size();k++) {
	double der_value;
	int t_plus = (t+1)%NT;
	int t_minus = (t-1+NT)%NT;
	double D_plus = D_List.distr_list[t_plus].distr[k];
	double D_0 = D_List.distr_list[t].distr[k];
	double D_minus = D_List.distr_list[t_minus].distr[k];

	if(mode==0) { //central derivative
	  der_value = (D_plus-D_minus)/2.0;
	}
	else if(mode==1) { //forward derivative
	  der_value = D_plus-D_0;
	}

	else if(mode==-1) { //backward derivative
	  der_value = D_0 - D_minus;
	}

	else crash("In static distr_t_list::derivative mode: "+to_string(mode)+" not yet implemented");

	RET.distr_list[t].distr.push_back(der_value);

      }


    }

    return RET;
  }
 



  bool UseJack;
  vector<distr_t> distr_list;
};

//operator overloading




  
distr_t_list operator*(const distr_t_list& A, const distr_t_list& B);

distr_t_list operator+(const distr_t_list& A,const distr_t_list& B);

distr_t_list operator-(const distr_t_list& A, const distr_t_list& B);

distr_t_list operator/(const distr_t_list& A,const distr_t_list& B) ;

Vfloat operator%(const distr_t_list& A,const distr_t_list& B);

////////////////////////////////

distr_t_list operator*(const distr_t_list& A, const distr_t& B);
distr_t_list operator*(const distr_t& B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list& A,const distr_t& B);
distr_t_list operator+(const distr_t& B, const distr_t_list& A);

distr_t_list operator-(const distr_t_list& A,const distr_t& B);

distr_t_list operator-(const distr_t& B,const distr_t_list& A);

distr_t_list operator/(const distr_t_list& A,const distr_t& B);

distr_t_list operator/(const distr_t& B,const distr_t_list& A);
  

Vfloat operator%(const distr_t_list& A,const distr_t& B);

Vfloat operator%(const distr_t& B, const distr_t_list& A);
///////////////////////////////


distr_t_list operator*(const distr_t_list& A, double B);
distr_t_list operator*(double B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list& A, double B);

distr_t_list operator+(double B, const distr_t_list& A);

distr_t_list operator-(const distr_t_list& A, double B);

distr_t_list operator-(double B, const distr_t_list& A);

distr_t_list operator/(const distr_t_list& A, double B);
  
distr_t_list operator/(double B, const distr_t_list& A);

distr_t_list operator+(const distr_t_list &A, const Vfloat &B);
distr_t_list operator-(const distr_t_list &A, const Vfloat& B);

distr_t_list operator+(const Vfloat &B, const distr_t_list &A);
distr_t_list operator-(const Vfloat& B, const distr_t_list& A);


distr_t_list operator*(const distr_t_list &A, const Vfloat& B);

distr_t_list operator*(const Vfloat& B, const distr_t_list& A);


//end operator overloading

class Jackknife {

 public:
  Jackknife() {}
  Jackknife(int block_size_max, int Njacks) {
    ThermalMeasToErase=0;
    this->block_size_max = block_size_max;
    this->Njacks= Njacks;
    Verbose_jack=0;
    ReturnBlockSizeMax=0;
    Enable_fractional_jackknife= true;
    }
  distr_t DoJack(int , ...);
  distr_t DoJack(function<double(const Vfloat&)> F, int Nobs, ... );

  


  int ThermalMeasToErase;
  int block_size_max;
  int Verbose_jack;
  int Njacks;
  bool ReturnBlockSizeMax;
  bool Enable_fractional_jackknife;


} ;




class Bootstrap {

 public:
  Bootstrap() {}
  Bootstrap(int Nboots, int seed, int size) : RN(seed, size-1) { 
    this->Nboots= Nboots;
    this->size=size;
    this->generate_data();
    
  }
  distr_t DoBoot(int Nobs, ...);
  distr_t DoBoot(function<double(const Vfloat&)> F, int Nobs, ... );
 
  
  int Nboots;
  int size;
  vector<vector<int>> bootstrap_data;

  private:
  RandomMersenne RN;
  void generate_data() {
    bootstrap_data.resize(Nboots);
    for(int iboot=0; iboot<Nboots;iboot++)
      for(int el=0; el < size; el ++) bootstrap_data[iboot].push_back( RN());
  }

} ;


distr_t Get_id_jack_distr(int N);
distr_t_list Get_id_jack_distr_list(int size, int N);
distr_t Get_id_distr(int N, bool UseJack);
distr_t_list Get_id_distr_list(int size, int N, bool UseJack);
distr_t AIC( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas, bool NEIL=false , const vector<double> &mult= vector<double>(1,0) );
distr_t AIC_lin( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas, bool NEIL=false, const vector<double> &mult= vector<double>(1,0));
distr_t BMA_uniform( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas, double ch2_th);
distr_t BMA_Eq_29( const vector<distr_t> &VAL);

distr_t_list EXP_DL(const distr_t_list &A);
distr_t_list EXPT_DL(const distr_t_list &A);
distr_t_list COSH_DL(const distr_t_list &A);
distr_t_list SINH_DL(const distr_t_list &A);
distr_t_list SQRT_DL(const distr_t_list &A);
distr_t_list LOG_DL(const distr_t_list &A);
distr_t_list POW_DL(const distr_t_list &A, int n);

distr_t_list EXPT_D(const distr_t &A, int T);
distr_t EXP_D(const distr_t &A);
distr_t COSH_D(const distr_t &A);
distr_t SINH_D(const distr_t &A);
distr_t SQRT_D(const distr_t &A);
distr_t LOG_D(const distr_t &A);
distr_t POW_D(const distr_t &A, int n);





#endif
