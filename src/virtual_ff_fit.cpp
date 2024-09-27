#include "../include/virtual_ff_fit.h"

using namespace std;

const int nsteps= 100;
const double step_size= 1.0/(double)nsteps;
double rk_guess_V = pow(493.677/775.4,2)  ; //Mk^2 /Mrho^2
double rq_guess_V = pow(493.677/892 , 2); //Mk^2/Mk*^2(892)
double rk_guess_A = pow(493.677/775.4,2) ; //Mk^2/Mrho^2
double rq_guess_H2 = 1.0;
double rq_guess_A = pow(493.677/1253 ,2); //Mk^2/Mk*^2(1270)
double phi2_mass_ratio = pow(493.677/1019,2);
const bool PURE_VMD=0;
double rk2_Mk2= 1.9540;
double fk_mk= 0.3140;

class Yff{

public:
  Yff() : Za_ov_Zv(1.0), Za_ov_Zv_err(0.0) {}
  double Mpi, fp, ff, xk, xq, ff_err, Za_ov_Zv, Za_ov_Zv_err;
  
};

class Xff_VMD{

public:
  Xff_VMD() {}
  Xff_VMD(const Vfloat &par) {
    if(par.size() != 5) {
      cout<<"In class Xff_VMD, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    a0=par[0];
    ampl=par[1];
    rk= par[2];
    rq=par[3];
    Za_ov_Zv=par[4];
  }

  double a0, ampl, rk,rq, Za_ov_Zv;
  


};

class Xff_ChPT{
  
public:
  Xff_ChPT() {}
  Xff_ChPT(const Vfloat &par) {
    if(par.size() != 6) {
      cout<<"In class Xff_ChPT, invalid call to constructor"<<endl;
      crash("par_size: +"+to_string(par.size()));
    }
    a0=par[0];
    ak= par[1];
    aq= par[2];
    akq= par[3];
    a2kq=par[4];
    Za_ov_Zv=par[5];
  }
  double a0, ak, aq, akq, a2kq, Za_ov_Zv;
};

 

void Fit_virtual_FF_VMD( vector<function<double(double, double)>> &fit_func,  distr_t &A_ret, distr_t &rk_ret, distr_t &rq_ret, const vector<distr_t> &FF, distr_t &f_p,distr_t &m_p, distr_t& Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag, string Meson, bool UseJack, bool ConstFit, int t_fit) {

  int Nmeas= FF.size();
  int njacks = f_p.size();
  if(f_p.size() != m_p.size()) crash(" Number of jackknives are not the same for fp and mp");
  bool verbose=1;


  cout<<"Nmeas: "<<Nmeas<<endl;
  cout<<"Njacks: "<<njacks<<endl;
  //perform bootstrap fit

  double L, m_l, T;
  Read_pars_from_ensemble_tag(Ens_tag, m_l, L, T);
  LatticeInfo L_info("LOCAL"); //CURRENT_TYPE == LOCAL


  bootstrap_fit<Xff_VMD,Yff> bf(njacks);


  bf.Set_number_of_measurements(FF.size());
  bf.Set_verbosity(verbose);
  bf.set_warmup_lev(1);
  

  //Add parameters
  double a0, a0_err;
  if(ff_type=="H1" || ff_type == "H2") {a0=0.2; a0_err=0.01;} 
  else if(ff_type=="FA")  {a0=0.03; a0_err=0.003;} 
  else if(ff_type=="FV") {a0=0.08*Za_ov_Zv.ave(); a0_err=0.002*Za_ov_Zv.ave();}
  else crash("string : "+ff_type+" does not name a form factor");
  bf.Add_par("a0",a0, a0_err);
  bf.Add_par("ampl",0.2  , 0.03);
  if(W=="A") {
    //bf.Add_prior_par("rk", rk_guess_A, rk_guess_A/10);   //COMMENT TO AVOID PRIOR ON RK RQ
    //bf.Add_prior_par("rq", rq_guess_A, rq_guess_A/10);   //COMMENT TO AVOID PRIOR ON RK RQ
    bf.Add_par("rk", rk_guess_A, rk_guess_A/10);
    if(ff_type != "H2") {
    bf.Add_par("rq", rq_guess_A, rq_guess_A/10);
    }
    else bf.Add_par("rq", rq_guess_A, rq_guess_A/10);
  
    if(ff_type=="H2") {  bf.Set_limits("rk", 0.0, 1.0);  bf.Set_limits("rq", -0.5, 2.0*rq_guess_A); }
    else  { bf.Set_limits("rk", 0.0, 3.0*rk_guess_A); bf.Set_limits("rq", -0.3, 6.0*rq_guess_A);} //4.0*rq_guess_A;
  }
  else if(W=="V") {
    
    //bf.Add_prior_par("rk", rk_guess_V, rk_guess_V/10);  //COMMENT TO AVOID PRIOR ON RK RQ
    //bf.Add_prior_par("rq", rq_guess_V, rq_guess_V/10);  //COMMENT TO AVOID PRIOR ON RK RQ
    bf.Add_par("rk", rk_guess_V, rk_guess_V/10);
    bf.Add_par("rq", rq_guess_V, rq_guess_V/10);
    bf.Set_limits("rk", 0.0, 3.0*rk_guess_V);
    bf.Set_limits("rq", 0.0, 3.0*rq_guess_V);
    
  }
  else crash("string W: "+W+" is not axial or vector");


  //fix n release
  //#######################################################################
  bf.Add_prior_par("Za_ov_Zv", Za_ov_Zv.ave(), Za_ov_Zv.err());
  if(W != "V") bf.Fix_par("Za_ov_Zv", 1.0);

  if(!PURE_VMD) {
    //bf.Fix_n_release("rk");        //COMMENT TO AVOID PRIOR ON RK RQ
    //bf.Fix_n_release("rq");        //COMMENT TO AVOID PRIOR ON RK RQ
  }
  if(W=="V") bf.Fix_n_release("Za_ov_Zv");


  if(PURE_VMD) {
    bf.Fix_par("a0",0.0);
    if(W=="V") {
      bf.Fix_par("rk", rk_guess_V);
      bf.Fix_par("rq", rq_guess_V);
    }
    else {
    bf.Fix_par("rk", rk_guess_A);
    bf.Fix_par("rq", rq_guess_A);
    }
  }

  //fix params to make test
  if(ff_type=="FV" || ff_type=="FA" ||  ff_type =="H1" || ff_type =="H2")   bf.Fix_par("a0",0.0);
  
 


  
  bf.ansatz =  [=](const Xff_VMD &p, const Yff &ip) -> double {
		 double val =  p.a0 + p.ampl/((1.0-pow(ip.xk,2)*p.rk)*(1.0- pow(ip.xq,2)*p.rq));
        
		 if(ff_type=="H2") val =p.a0 + p.ampl*(1.0/((1.0-pow(ip.xk,2)*p.rk)*(1.0-pow(ip.xq,2)*p.rq)))     ;
		 return val;
  };
  bf.measurement = [=](const Xff_VMD& p,const Yff& ip) -> double {
    return p.Za_ov_Zv*ip.ff;
  };
  bf.error =  [=](const Xff_VMD& p,const  Yff &ip) -> double {
    return p.Za_ov_Zv*ip.ff_err;
  };


  //where to store fitted parameters
  boot_fit_data<Xff_VMD> Bt_fit;

  //input parameters
  vector<vector<Yff>> data(njacks);
    for(auto &data_iboot: data) data_iboot.resize(Nmeas);
  
    for(int ijack=0;ijack<njacks;ijack++) {
      for(int imeas=0;imeas<Nmeas;imeas++) {
        data[ijack][imeas].Mpi = m_p.distr[ijack];
	data[ijack][imeas].fp = f_p.distr[ijack];
	data[ijack][imeas].ff = FF[imeas].distr[ijack];
	data[ijack][imeas].ff_err = FF[imeas].err();
	data[ijack][imeas].xk = fabs(mom[imeas].virt()/m_p.distr[ijack]) ;
	data[ijack][imeas].xq = sqrt( 1.0 + pow(data[ijack][imeas].xk,2) - mom[imeas].x_gamma(m_p).distr[ijack]) ;
      } 
    }
   
 

  
  //add prior
      for(int ijack=0;ijack<njacks;ijack++) {
	bf.ib= &ijack;
        //bf.Append_to_prior("rk", rk_guess_A, rk_guess_A );            //COMMENT TO AVOID PRIOR ON RK RQ
	if(W=="A")  {
	  //bf.Append_to_prior("rq", rq_guess_A, rq_guess_A );          //COMMENT TO AVOID PRIOR ON RK RQ
	  bf.Append_to_prior("Za_ov_Zv", 1.0, 1.0);
	}
	else if(W=="V") {
	  //bf.Append_to_prior("rq", rq_guess_V, rq_guess_V );          //COMMENT TO AVOID PRIOR ON RK RQ
	  bf.Append_to_prior("Za_ov_Zv", Za_ov_Zv.distr[ijack], Za_ov_Zv.err());
	}
	else crash("string W: "+W+" is not axial nor vector");
      }
    


  //fit
  bf.Append_to_input_par(data);
  Bt_fit =bf.Perform_bootstrap_fit();


  //print fitted params
  Vfloat a0_boot, ampl_boot, rk_boot, rq_boot;
  for (auto &p: Bt_fit.par) {
    a0_boot.push_back(p.a0);
    ampl_boot.push_back(p.ampl);
    rk_boot.push_back(p.rk);
    rq_boot.push_back(p.rq);
  }

  cout<<"########PRINTING FIT PARAMS VMD FIT OF "<<ff_type<<"###########"<<endl;
  cout<<"A: "<<Boot_ave(a0_boot)<<" +- "<<Boot_err(a0_boot,UseJack)<<endl;
  cout<<"B: "<<Boot_ave(ampl_boot)<<" +- "<<Boot_err(ampl_boot,UseJack)<<endl;
  cout<<"C: "<<Boot_ave(rk_boot)<<" +- "<<Boot_err(rk_boot,UseJack)<<endl;
  cout<<"D: "<<Boot_ave(rq_boot)<<" +- "<<Boot_err(rq_boot,UseJack)<<endl;
  cout<<"##########FIT PARAMS PRINTED###########"<<endl;

  //params printed


  //define lambda functions and return fit parameters
  for(int ijack=0;ijack<njacks;ijack++) {
    auto F = [=](double xk, double xq) -> double {
	       double val= Bt_fit.par[ijack].a0 + Bt_fit.par[ijack].ampl/((1.0-pow(xk,2)*Bt_fit.par[ijack].rk)*(1.0 -pow(xq,2)*Bt_fit.par[ijack].rq));
	       if(ff_type=="H2") val = Bt_fit.par[ijack].a0+  Bt_fit.par[ijack].ampl*(1.0/((1.0 -pow(xk,2)*Bt_fit.par[ijack].rk)*(1.0 -pow(xq,2)*Bt_fit.par[ijack].rq)));
	       return val;
	     };
    
    fit_func.push_back(F);
    
    //ret parameters
    A_ret.distr.push_back( Bt_fit.par[ijack].ampl);
    rk_ret.distr.push_back( Bt_fit.par[ijack].rk);
    rq_ret.distr.push_back( Bt_fit.par[ijack].rq);
  }


  
  //compute form factors from ChPT
  function<double(double,double)> ChPT_H1,ChPT_H2,ChPT_FA,ChPT_FV,F_to_Print;
  Compute_ChPT_form_factors(ChPT_H1,ChPT_H2, ChPT_FA, ChPT_FV);
  if(ff_type=="FA") F_to_Print= ChPT_FA;
  else if(ff_type=="FV") F_to_Print= ChPT_FV;
  else if(ff_type=="H1") F_to_Print=ChPT_H1;
  else if(ff_type=="H2") F_to_Print=ChPT_H2;
  else crash("cannot interpret form factor type: "+ff_type);
  
  //print data
  //set a grid in the xk, xq plane

  Vfloat xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list;
  for(int dxk=0; dxk < nsteps; dxk++) {
    for(int dxq=0; dxq < nsteps; dxq++) {
      double xxk = step_size*dxk;
      double xxq = step_size*dxq;
      distr_t func(UseJack);
      for(int ijack=0;ijack<njacks;ijack++) {
	func.distr.push_back( fit_func[ijack](xxk,xxq));
      }
      xk_list.push_back(xxk);
      xq_list.push_back(xxq);
      FF_list.push_back(func.ave());
      FF_err_list.push_back(func.err());
      FF_ChPT_list.push_back(F_to_Print(xxk,xxq));
    }
  }


 
  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag);
  if(ConstFit) {
    Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag+"/"+ff_type+".dat", "", "#xk   xq  val    err   ChPT pred" );
  }
  else  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_VMD_"+Ens_tag+"/"+ff_type+"time_"+to_string(t_fit)+".dat", "", "#xk   xq  val    err    ChPT_pred" );
  
  
 

  return;

}







void Fit_virtual_FF_ChPT(vector<function<double(double, double)>> &fit_func,  distr_t &a0_ret, distr_t &ak_ret, distr_t &aq_ret, const vector<distr_t> &FF, distr_t &f_p, distr_t &m_p, distr_t &Za_ov_Zv, vector<pt3_momenta> &mom, string ff_type, string W, string Ens_tag, string Meson, bool UseJack, bool ConstFit, int t_fit) {

  int Nmeas= FF.size();
  int njacks = f_p.size();
  if((f_p.size() != m_p.size()) || (f_p.size() != Za_ov_Zv.size())) crash(" Number of jackknives are not the same for fp,mp and RC ZA/ZV");
  bool verbose=1;

 cout<<"Nmeas: "<<Nmeas<<endl;
 cout<<"Njacks: "<<njacks<<endl;
  //perform bootstrap fit

 
  bootstrap_fit<Xff_ChPT,Yff> bf(njacks);


  bf.Set_number_of_measurements(FF.size());
  bf.Set_verbosity(verbose);
  

  //Add parameters
  if(ff_type=="H1")  { bf.Add_par("a0", 0.2, 0.004); bf.Add_par("ak", 0.12 , 0.008); bf.Add_par("aq", 0.09 ,0.001);}
  else if(ff_type == "H2") { bf.Add_par("a0", 0.2, 0.005); bf.Add_par("ak",0.23 ,0.003); bf.Add_par("aq",-0.1 ,0.002);}
  else if(ff_type=="FA") {bf.Add_par("a0",0.032, 0.002); bf.Add_par("ak",0.036 ,0.001); bf.Add_par("aq",-0.0002 ,0.00001);}
  else if(ff_type=="FV") {bf.Add_par("a0",0.074*Za_ov_Zv.ave(), 0.001*Za_ov_Zv.ave()); bf.Add_par("ak",0.045*Za_ov_Zv.ave() ,0.0001*Za_ov_Zv.ave()); bf.Add_par("aq",0.02*Za_ov_Zv.ave() ,0.001*Za_ov_Zv.ave());}
  else crash("string : "+ff_type+" does not name a form factor");
  bf.Add_par("akq", 0.1, 0.01);
  bf.Add_par("a2kq",0.3,0.01);
  bf.Fix_par("akq",0.0);
  //bf.Fix_par("akq", 0.0);
  //bf.Fix_par("aq",0.0);
  //bf.Fix_par("akq",0.0);
  //bf.Fix_par("ak",0.0);
  bf.Fix_par("a2kq",0.0);                                    //                    YOU FIX HERE TO ZERO THE PARAMETER a2kq

 
  //fix n release Za/Zv
  bf.Add_prior_par("Za_ov_Zv", Za_ov_Zv.ave(), Za_ov_Zv.err());
  if(W != "V") bf.Fix_par("Za_ov_Zv", 1.0);
  if(W =="V") bf.Fix_n_release("Za_ov_Zv");

  

  // if(ff_type=="H1" || ff_type=="H2") {bf.Fix_par("ak",0.0);bf.Fix_par("aq",0.0); bf.Fix_par("akq",0.0);}
  bf.ansatz =  [=](const Xff_ChPT &p, const Yff &ip) -> double {
    return  p.a0 + p.ak*pow(ip.xk,2) + p.aq*pow(ip.xq,2)+ p.akq*ip.xq*ip.xk +  p.a2kq*pow(ip.xq*ip.xk,2)  ;
  };
  bf.measurement = [=](const Xff_ChPT& p,const Yff& ip) -> double {
    return p.Za_ov_Zv*ip.ff;
  };
  bf.error =  [=](const Xff_ChPT& p,const  Yff &ip) -> double {
    return p.Za_ov_Zv*ip.ff_err;
  };


  //where to store fitted parameters
  boot_fit_data<Xff_ChPT> Bt_fit;

  //input parameters
  vector<vector<Yff>> data(njacks);
  for(auto &data_iboot: data) data_iboot.resize(Nmeas);

  
    for(int ijack=0;ijack<njacks;ijack++) {
      for(int imeas=0;imeas<Nmeas;imeas++) {
        data[ijack][imeas].Mpi = m_p.distr[ijack];
	data[ijack][imeas].fp = f_p.distr[ijack];
	data[ijack][imeas].ff = FF[imeas].distr[ijack];
	data[ijack][imeas].ff_err = FF[imeas].err();
	data[ijack][imeas].xk = fabs(mom[imeas].virt()/m_p.distr[ijack]) ;
	data[ijack][imeas].xq = sqrt( 1.0 + pow(data[ijack][imeas].xk,2) - mom[imeas].x_gamma(m_p).distr[ijack]) ;
      } 
    }


     //add prior 
  for(int ijack=0;ijack<njacks;ijack++) {
     bf.ib= &ijack;
     if(W=="A")  {  bf.Append_to_prior("Za_ov_Zv", 1.0, 1.0);}
     else if(W=="V") { bf.Append_to_prior("Za_ov_Zv", Za_ov_Zv.distr[ijack], Za_ov_Zv.err());}
     else crash("string W: "+W+" is not axial nor vector");
  }

    
  //fit
  bf.Append_to_input_par(data);
  Bt_fit =bf.Perform_bootstrap_fit();


   //print fitted params
  Vfloat a0_boot, ak_boot, aq_boot, akq_boot, a2kq_boot;
  for (auto &p: Bt_fit.par) {
    a0_boot.push_back(p.a0);
    ak_boot.push_back(p.ak);
    aq_boot.push_back(p.aq);
    akq_boot.push_back(p.akq);
    a2kq_boot.push_back(p.a2kq);
  }

  cout<<"########PRINTING FIT PARAMS CHPT FIT OF "<<ff_type<<"###########"<<endl;
  cout<<"a0: "<<Boot_ave(a0_boot)<<" +- "<<Boot_err(a0_boot,UseJack)<<endl;
  cout<<"ak: "<<Boot_ave(ak_boot)<<" +- "<<Boot_err(ak_boot,UseJack)<<endl;
  cout<<"aq: "<<Boot_ave(aq_boot)<<" +- "<<Boot_err(aq_boot,UseJack)<<endl;
  cout<<"akq: "<<Boot_ave(akq_boot)<<" +- "<<Boot_err(akq_boot,UseJack)<<endl;
  cout<<"a2kq: "<<Boot_ave(a2kq_boot)<<" +- "<<Boot_err(a2kq_boot,UseJack)<<endl;
  cout<<"##########FIT PARAMS PRINTED###########"<<endl;

  //params printed


  //define lambda functions and return parameters
  for(int ijack=0;ijack<njacks;ijack++) {
    auto F = [=](double xk, double xq) -> double { return Bt_fit.par[ijack].a0 + Bt_fit.par[ijack].ak*pow(xk,2) + Bt_fit.par[ijack].aq*pow(xq,2) +Bt_fit.par[ijack].akq*xk*xq + Bt_fit.par[ijack].a2kq*pow(xk*xq,2);};
    fit_func.push_back(F);

    //ret parameters
    a0_ret.distr.push_back( Bt_fit.par[ijack].a0);
    ak_ret.distr.push_back( Bt_fit.par[ijack].ak);
    aq_ret.distr.push_back( Bt_fit.par[ijack].aq);
  }


  //compute form factors from ChPT
  function<double(double,double)> ChPT_H1,ChPT_H2,ChPT_FA,ChPT_FV,F_to_Print;
  Compute_ChPT_form_factors(ChPT_H1,ChPT_H2, ChPT_FA, ChPT_FV);
  if(ff_type=="FA") F_to_Print= ChPT_FA;
  else if(ff_type=="FV") F_to_Print= ChPT_FV;
  else if(ff_type=="H1") F_to_Print=ChPT_H1;
  else if(ff_type=="H2") F_to_Print=ChPT_H2;
  else crash("cannot interpret form factor type: "+ff_type);
  
  //print data
  //set a grid in the xk, xq plane
 

  Vfloat xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list;
  for(int dxk=0; dxk < nsteps; dxk++) {
    for(int dxq=0; dxq < nsteps; dxq++) {
      double xxk = step_size*dxk;
      double xxq = step_size*dxq;
      distr_t func(UseJack);
      for(int ijack=0;ijack<njacks;ijack++) {
	func.distr.push_back( fit_func[ijack](xxk,xxq));
      }
      xk_list.push_back(xxk);
      xq_list.push_back(xxq);
      FF_list.push_back(func.ave());
      FF_err_list.push_back(func.err());
      FF_ChPT_list.push_back(F_to_Print(xxk,xxq));
    }
  }


  boost::filesystem::create_directory("../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag);
   
  if(ConstFit) {
    Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag+"/"+ff_type+".dat", "", "#xk   xq  val    err     ChPT_pred" );
  }
  else  Print_To_File( {}, {xk_list, xq_list, FF_list, FF_err_list, FF_ChPT_list}, "../data/form_factors/"+Meson+"/virtual_FF_fit_ChPT_"+Ens_tag+"/"+ff_type+"time_"+to_string(t_fit)+".dat", "", "#xk   xq  val    err     ChPT_pred" );
  
 

  return;

}
