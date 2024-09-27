#ifndef __bootstrap_fit__
#define __bootstrap_fit__


#include "numerics.h"
#include "stat.h"

using namespace std;

template <class T>
class boot_fit_data {

public:
  
  boot_fit_data() {}
  void ch2_ave() {
    double ch2=0.0, err_ch2=0.0;
    double N= (double)chi2.size();
    for(auto &c :chi2) {ch2+=c/N; err_ch2 += c*c/N;}
    cout<<"average bootstrap chi2: "<<ch2<<" +- "<<sqrt( ((N-1.0)/N)*(err_ch2- ch2*ch2))<<endl;
    return;
  }
  double get_ch2_ave() {
    double ch2=0.0;
    double N= (double)chi2.size();
    for(auto &c :chi2) ch2+=c/N;
    return ch2;
  }

  double get_ch2_err() {
    distr_t ch2_distr;
    for(int i=0;i<(signed)chi2.size();i++) ch2_distr.distr.push_back( chi2[i]);
    return ch2_distr.err();

  }
  Vfloat chi2;
  Vfloat EDM;
  Vfloat N_it;
  vector<bool> IsValid;
  vector<T> par;

};





template <class T1, class T2> 

class bootstrap_fit : public ROOT::Minuit2::FCNBase  {
  
public:
  bootstrap_fit():  theErrorDef(1.0), NumberOfMeasurements(0), warm_up(0) {
    this->PATH="chi2.out"; this->verbose=0; this->Use_Cov_Matrix=false;
  }
  bootstrap_fit(int nb) :  Input_pars(nb), theErrorDef(1.0), NumberOfMeasurements(0), nboots(nb), warm_up(0) { this->PATH="chi2.out"; this->verbose=0; this->Use_Cov_Matrix=false;}
 
  virtual ~bootstrap_fit() {}
  function<double(const T1& p, const T2 &ip)> ansatz;
  function<double(const T1& p, const T2 &ip)> measurement;
  function<double(const T1& p, const T2 &ip)> error;
  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const Vfloat& par) const;

  boot_fit_data<T1> Perform_bootstrap_fit();
  
  
  void setErrorDef(double def) {theErrorDef = def;}
  void Add_par(string Name, double val, double err);
  void Set_par_val(string Name, double val);
  void Set_par_val(string Name, double val, double err);
  void Add_prior_par(string Name, double val, double err);
  void Add_prior_pars(const vector<string>& Names);
  void Append_to_prior(string Name, double val, double err);
  void Fix_par(string Name, double val);
  void Release_par(string Name);
  void Fix_par(string Name) { Fix_par(Name, nan("1"));};
  void Append_to_input_par(vector<vector<T2>> values) {this->Input_pars = values;}
  void Append_to_input_par(vector<T2> val) { this->Input_pars[*ib] = val;}
  void Append_to_input_par(T2 val);
  void Clear_input_pars() {Input_pars.clear(); Input_pars.resize(nboots);}
  void Clear_prior(string A);
  void Clear_priors() { for (auto &pr: Priors) pr.second.clear();}
  void Delete_priors() {Priors.clear();}
  void Set_number_of_measurements(int Nens) { this->NumberOfMeasurements = Nens;}
  int Get_number_of_measurements() {return this->NumberOfMeasurements;}
  void Set_verbosity(int a) { this->verbose=a;}
  void Set_nboots(int nboots) { this->nboots=nboots;}
  void Fix_n_release(string A, double val) {this->To_release.insert(make_pair(A, val));}
  void Fix_n_release(string A) {Fix_n_release(A, nan("1"));}
  void Fix_n_release(const vector<string> &A) { for( auto &p:A) Fix_n_release(p);}
  int Get_number_of_fit_pars() { return PNames.size()  - Fixed_pars.size() -To_release.size();}
  void Set_limits(string Name, double val_1, double val_2) {

    PList.SetLimits(Name, val_1,val_2);
    return;  
  }
  int P(string A) const ;
  void set_warmup() {warm_up=true;};
  void set_warmup_lev(int ilev) {if(ilev <0) crash("cannot set warm up lev to negative number");warm_up=ilev;}
  void Set_print_path(string A) { this->PATH = A;}
  void Add_covariance_matrix(Eigen::MatrixXd M) {
    this->Use_Cov_Matrix=true;
    if( M.rows() != NumberOfMeasurements) crash("Size of covariant matrix and number of measurements do not match");
    this->Cov_matrix_inv= M.inverse();
  }
  void Disable_correlated_fit() { this->Use_Cov_Matrix=false;}


  
 
  int* ib;

 private:
  ROOT::Minuit2::MnUserParameters PList; 
  map<string, double> To_release;
  map<string, int> PNames;
  map<string, double> Fixed_pars;
  vector<vector<T2>> Input_pars;
  map<string, VPfloat> Priors;
  double theErrorDef;
  int NumberOfMeasurements;
  int nboots;
  int verbose;
  bool PRINT;
  string TAG;
  string PATH;
  int warm_up;
  bool Use_Cov_Matrix;
  Eigen::MatrixXd Cov_matrix_inv;
 
  
};


template <class T1, class T2> 
double bootstrap_fit<T1, T2>::operator()(const Vfloat& par) const {
    int pos=0;
    double chi2=0.0;
    ofstream PrintChi;
    if(PRINT) {
      PrintChi.open(PATH, ofstream::app); PrintChi.precision(10);
      PrintChi<<"#############################################"<<endl;
      PrintChi<<this->TAG<<endl;
    }


    if(Use_Cov_Matrix) {
      for(int i=0;i<NumberOfMeasurements;i++) {
	T1 pp(par);
	double ansatz_i= this->ansatz(pp, Input_pars[*ib][i]);
	double measurement_i = this->measurement(pp, Input_pars[*ib][i]);
	double res_i = 0;
	for(int j=0;j<NumberOfMeasurements;j++) {
	  double ansatz_j = this->ansatz(pp, Input_pars[*ib][j]);
	  double measurement_j = this->measurement(pp, Input_pars[*ib][j]);
	  double res= (ansatz_i -measurement_i)*Cov_matrix_inv(i,j)*(ansatz_j-measurement_j);
	  chi2+=res;
	  res_i+=res;
	}
	if(PRINT) PrintChi<<i<<setw(20)<<ansatz_i<<setw(20)<<measurement_i<<setw(20)<<res_i<<endl;
      }
    }
    
    else {
      while(pos<NumberOfMeasurements)  {
      T1 pp(par);
      // for(unsigned int pi=0;pi<par.size();pi++) if(isnan(par[pi])) crash("par nr: "+to_string(pi)+" is nan, bootstrap: "+to_string(*ib));
      double ansatz = this->ansatz(pp, Input_pars[*ib][pos]);
      double measurement = this->measurement(pp, Input_pars[*ib][pos]);
      double error = this->error(pp, Input_pars[*ib][pos]);
          
      double res = pow( (ansatz-measurement)/error,2);
      if(PRINT) PrintChi<<pos<<setw(20)<<ansatz<<setw(20)<<measurement<<setw(20)<<error<<setw(20)<<res<<endl;
      chi2+=res;
      pos++;
      }
    }

    if(PRINT) PrintChi<<"chi2 w.o. priors: "<<chi2<<endl;
    //Add Gaussian Prior
    for(auto const& [key, val] : this->Priors) {
      double ch2_pr = pow((par[P(key)]- val[*ib].first)/val[*ib].second,2);
      if(PRINT) PrintChi<<key<<setw(20)<<par[P(key)]<<setw(20)<<val[*ib].first<<setw(20)<<val[*ib].second<<setw(20)<<ch2_pr<<endl;
      chi2 += ch2_pr;
    }
    if(PRINT) PrintChi<<"chi_2: "<<chi2<<endl;
    if(PRINT) PrintChi.close();
    return chi2;
}



template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Append_to_input_par(T2 val) {

  if(Input_pars[*ib].size() >= this->NumberOfMeasurements) crash("In bootstrap_fit::Append_to_input_par, data size larger than NumberOfMeasurements");
  else Input_pars[*ib].push_back(val);
  return;
}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Clear_prior(string A) {
 

  for(auto & [key,val] : Priors) {
    if(key == A) val.clear();
  }
  return;
}

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Append_to_prior(string Name, double val, double err) {

  map<string, VPfloat>::iterator it;
  it = Priors.find(Name);
  if(it == Priors.end()) crash("Could not find the prior parameter "+Name+" in Prior map");
  else {
    it->second.push_back(make_pair(val,err));
  }
  return;
  

}



template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_par(string Name, double val, double err) {
    PList.Add(Name, val, err);
    int Position= PNames.size();
    PNames.insert( pair<string,int>(Name, Position));
    return; 
  }

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Set_par_val(string Name, double val) {
    PList.SetValue(Name, val);
    return; 
  }

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Set_par_val(string Name, double val, double err) {
    PList.SetValue(Name, val);
    PList.SetError(Name,val);
    return; 
  }


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_prior_pars(const vector<string> &Names) {

  
    for(unsigned int i=0; i<Names.size(); i++) {
      this->PList.Add(Names[i], 0.0, 1.0);
      int Position=PNames.size();
      PNames.insert( pair<string, int>(Names[i], Position));
      VPfloat start;
      Priors.insert(pair<string, VPfloat>(Names[i], start));
    }
    return;
}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Add_prior_par(string Name, double val, double err) {

   
  PList.Add(Name, val, err);
  int Position=PNames.size();
  PNames.insert( pair<string, int>(Name, Position));
  VPfloat start;
  Priors.insert(pair<string, VPfloat>(Name, start));
  
  return;
}

template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Fix_par(string Name, double val) {
  map<string,double>::iterator it;
  it = Fixed_pars.find(Name);
  if(it != Fixed_pars.end()) {
    it->second = val;
    return;
  }
  
  Fixed_pars.insert(make_pair(Name, val)); 
  return;  
}


template <class T1, class T2> 
void bootstrap_fit<T1,T2>::Release_par(string Name) {

  Fixed_pars.erase(Name);
  return;  
}


template <class T1, class T2> 
boot_fit_data<T1> bootstrap_fit<T1,T2>::Perform_bootstrap_fit() {

  boot_fit_data<T1> boot_result;

  int boot_index=0;
  ib = &boot_index;

  ROOT::Minuit2::MnMigrad migrad(*this, PList, 3);
  
  
  
  while(boot_index<nboots) {

     
    //fix parameters that must be fixed initially
    for(auto &name: To_release) {
      if(!isnan(name.second)) migrad.SetValue(name.first.c_str(), name.second);
      map<string, VPfloat>::iterator it;
      it= Priors.find(name.first);
      if( it != Priors.end()) migrad.SetValue(name.first.c_str(), it->second[*ib].first);
      migrad.Fix(name.first.c_str());
    }
    for(auto &name: Fixed_pars) {
      if(!isnan(name.second)) migrad.SetValue(name.first.c_str(), name.second);
      migrad.Fix(name.first.c_str());
    }
    

    
    //fit
    ROOT::Minuit2::FunctionMinimum chi2 = migrad();
    //warmup
    for(int i=0;i<warm_up;i++) chi2=migrad();

    if(verbose) { //print chi^2 details
      this->PRINT= true;
      this->TAG= "Fixed priors";
      Vfloat PP(PNames.size());
      for(auto const& [key,val]: PNames) PP[val] = chi2.UserState().Value(key);
      this->operator()(PP);
      this->PRINT = false;
    }



    //release them
    for(auto & [key, val] : To_release) migrad.Release(key.c_str());

    
    //refit with warm_up warmups
    for(int i=0;i< warm_up+1;i++)  chi2 = migrad();
    boot_result.chi2.push_back(chi2.Fval());
    boot_result.EDM.push_back(chi2.Edm());
    boot_result.N_it.push_back(chi2.NFcn());
    boot_result.IsValid.push_back(chi2.IsValid());
    Vfloat PP(PNames.size());
    for(auto const& [key,val]: PNames) PP[val] = chi2.UserState().Value(key);
    T1 P(PP);
    boot_result.par.push_back(P);

    if(verbose) { //print chi^2 details
      this->PRINT= true;
      this->TAG= "Priors released";
      this->operator()(PP);
      this->PRINT = false;
    }
    

    
    

    cout<<"Chi2 for bootstrap: "<<boot_index<<" = "<<chi2.Fval()<<endl;
    cout<<"Parameters value:"<<endl;
    for(unsigned int par_id=0; par_id<PNames.size();par_id++) {
      for(auto const& [key,val] : PNames) {
	if(val==(signed)par_id) cout<<key<<"\t"<<chi2.UserState().Value(key)<<endl;
      }
    }
    if(verbose) {
      cout<<"Exit status: "<<chi2.IsValid()<<endl;
      cout<<"EDM :"<<chi2.Edm()<<endl;
      cout<<"N_iterations :"<<chi2.NFcn()<<endl;
    }
    cout<<"#################"<<endl;


      


    boot_index++;
  }

  
  return boot_result;
}

template <class T1, class T2> 
int bootstrap_fit<T1,T2>::P(string A) const {

  int occurrences=0;
  int pos=0;

  for(auto const& [key, val] : this->PNames)
    if(key== A) { pos = val; occurrences++;}

  if(occurrences!= 1) crash("In bootstrap_fit<T1,T2>::P(string A), either par_name: "+A+" does not exist or it's present multiple times");

  return pos;
  
}
 





double Boot_ave(Vfloat& A);
Pfloat Boot_ave_err(Vfloat& A);
Pfloat Boot_ave_err(Vfloat& A, bool mode);
double Boot_err(Vfloat& A) ;
double Boot_err(Vfloat& A, bool mode);
Pfloat Boot_ave_err(VVfloat& A) ;
Pfloat Boot_ave_err(VVfloat& A, bool mode);
Pfloat Boot_ave_err(VVfloat& A, double resc, bool mode);
double Boot_ave(VVfloat& A) ;
double Boot_err(VVfloat& A);
double Boot_err(VVfloat& A, bool mode);
double Boot_err(VVfloat& A, double resc, bool mode);




































#endif

