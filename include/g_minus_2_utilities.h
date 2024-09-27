#ifndef  __g_minus_2_utilities__
#define  __g_minus_2_utilities__


#include "numerics.h"
#include "random.h"
#include "Corr_analysis.h"
#include "stat.h"
#include "LatInfo.h"
#include "T_min.h"
#include "input.h"
#include "highPrec.h"
#include "Bootstrap_fit.h"




using namespace std;



double kernel_K(double t,double MV);
double kernel_K_old(double t, double MV);
double der_kernel_K_W_win(double t, double MV);
double der_kernel_K_SD_win(double t, double MV);
double der_kernel_K(double t, double MV);
void Plot_kernel_K(int Npoints);
double Kernel_Pi_q2(double t, double Q, double a);
void Plot_Energy_windows_K();
double Zeta_function_laplacian_Luscher(double z);
void Zeta_function_zeroes(int Nzeros, Vfloat &res);
double tan_phi(double z);
double phi(double z);
double tan_phi_der(double z);
double phi_der(double z);
double phi_der_for_back(double z, int mode);
void Generate_free_corr_data();
void Get_spec_dens_free(const Vfloat &ams, string out_path);
void Compute_free_spectral_density(int Nc, double am, int reg, double step_size_erg, string dir_out);
void Compute_SD_window_Free();
void Compute_free_corr(double am, int Tmax);
double free_vector_corr_cont(int Nc, double am, double t);

class LL_functions{

public:

  LL_functions(VVfloat &phiD, VVfloat &phi_derD, Vfloat &sx_der, Vfloat &dx_der, Vfloat &sx_int, Vfloat &Dz, int Nresonances, Vfloat &L_zeroes) : Nres(Nresonances), sx_intervals(sx_int), Luscher_zeroes(L_zeroes) {
    //sanity checks
    if(phiD.size() != phi_derD.size() || phiD.size() != sx_der.size() || phiD.size() != dx_der.size() || phiD.size() != sx_int.size() || phiD.size() != Dz.size() ) crash("In constructor of LL_function size of vectors is not constant");

    int N= (signed)phiD.size();

    assert(N>0);

    for(int i=0;i<N;i++) {
      phi_spline_list.emplace_back(phiD[i].begin(), phiD[i].end(), sx_intervals[i], Dz[i],  sx_der[i], dx_der[i]);
      phi_der_spline_list.emplace_back(phi_derD[i].begin(), phi_derD[i].end(), sx_intervals[i], Dz[i]);
    }
    

    phi_spline= [&](double z) -> double {

		  bool find_pos=false;
		  double z2= z*z;
		  int pos=0;
		  for(unsigned int i=0; i<sx_intervals.size()-1;i++) {
		    if(sx_intervals[i+1] >= z2 && z2 >= sx_intervals[i] ) {pos=i;find_pos=true; break;}
		  }
		  if(!find_pos) pos= sx_intervals.size() -1;
		  return phi_spline_list[pos](z2);
		};


    phi_der_spline= [&](double z) -> double {

		  bool find_pos=false;
		  double z2= z*z;
		  int pos=0;
		  for(unsigned int i=0; i<sx_intervals.size()-1;i++) {
		    if(sx_intervals[i+1] >= z2 && z2 >= sx_intervals[i] ) {pos=i;find_pos=true; break;}
		  }
		  if(!find_pos) pos= sx_intervals.size() -1;
		  return phi_der_spline_list[pos](z2);
		};
    
  


  }


  double Vdual(double t, double m_rho, double Edual, double Rdual);
  double Gamma_rpp(double omega, double g_rho_pipi, double Mpi);
  double Gamma_rpp_der(double omega, double g_rho_pipi, double Mpi);
  double h(double omega, double g_rho_pipi, double Mpi);
  double h_prime(double omega, double g_rho_pipi, double Mpi);
  double h_second(double omega, double r_rho_pipi, double Mpi);
  double A_pipi_0(double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double Re_A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double Im_A_pipi(double omega, double g_rho_pipi, double Mpi);
  Pfloat A_pipi(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double F_pi_GS_mod(double omega, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double cot_d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double d_11(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double cot_d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double d_11_der(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double d_11_der_num(double k, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  double Amplitude(double k, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa);
  void Find_pipi_energy_lev(double L, double m_rho, double g_rho_pipi, double Mpi, double kappa,  Vfloat &res);
  double V_pipi(double t, double L, double m_rho, double g_rho_pipi, double Mpi, double kappa, Vfloat &Knpp);
  double V_pipi_infL(double t, double m_rho_infL, double g_rho_pipi_infL, double Mpi_infL, double kappa);
  void MLLGS_fit_to_corr(const distr_t_list &Corr,const distr_t &Mpi,const distr_t &a_distr, double L, distr_t &Edual, distr_t &Rdual, distr_t &Mrho, distr_t &grpp, int tmin, int tmax, string Tag);
  int Nres;
  Vfloat sx_intervals;
  Vfloat Luscher_zeroes;
  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> phi_spline_list;
  vector<boost::math::interpolators::cardinal_cubic_b_spline<double>> phi_der_spline_list;
  function<double(double)> phi_spline, phi_der_spline;

  

};








#endif
