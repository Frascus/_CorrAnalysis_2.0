#include "../include/T_min.h"


using namespace std;

double T_fit::operator()(const Vfloat& par) const {

   double ch2=0.0;
   for(int i_meas=0; i_meas < NumberOfMeasurements;i_meas++) {
     ch2+= pow( (this->ansatz(par,X[i_meas],i_meas)-Y[i_meas])/Y_err[i_meas],2);
   }

   return ch2;
}




fit_t_res T_fit::fit() {

  fit_t_res fit_result;

ROOT::Minuit2::MnUserParameters Param_List;
//add params
 int count_par=0;
 for(auto &par : init_pars) {
   string par_name="par"+to_string(count_par);
   if( (signed)init_par_errs.size() > count_par) 
   Param_List.Add(par_name , par, init_par_errs[count_par] );
   else Param_List.Add(par_name, par, fabs(par)/10.0);
   count_par++;
 }


 ROOT::Minuit2::MnMigrad migrad(*this, Param_List, 2);



ROOT::Minuit2::FunctionMinimum chi2 = migrad();

 
   
   fit_result.chi2 = chi2.Fval();
   fit_result.Iterations = chi2.NFcn();
   fit_result.Status = chi2.IsValid();
   for(int i=0;i<count_par;i++) {
     string retrieve_par_name="par"+to_string(i);
     fit_result.pars.push_back(chi2.UserState().Value(retrieve_par_name));
     //cout<<fit_result.pars[i]<<endl;
   }
   

   return fit_result;

}
