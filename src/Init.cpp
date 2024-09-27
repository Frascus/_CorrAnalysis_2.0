#include "../include/init.h"
#include "RC_WI_analysis.h"
#include "RC_analysis.h"
#include "chi_mag.h"


using namespace std;






MasterClass_analysis::MasterClass_analysis(string Path) {  //default constructor


  ifstream ReadInput(Path);

  //Read Header
  string Header;
  getline(ReadInput,Header);
  if(Header != "# Input file for the analysis of QCD+QED observables") crash("Input file does not have the expected header");

  string ReadLine;
  do {getline(ReadInput,ReadLine);} while(ReadLine=="" || ReadInput.eof());
  if(ReadLine != "#Analysis_Mode") crash("Error in input file. Expected: #Analysis_Mode, obtained: "+ReadLine+"|");
  ReadInput >> Analysis_Mode;


  
  if(Analysis_Mode=="Meson_masses") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }

  else if(Analysis_Mode=="Meson_masses_ov_X") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }

   
  else if(Analysis_Mode=="Meson_masses_twisted") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }

  else if(Analysis_Mode=="Meson_masses_twisted_adim") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }

  else if(Analysis_Mode=="Meson_masses_twisted_ov_X") {
    do {getline(ReadInput, ReadLine);} while(ReadLine=="" || ReadInput.eof());
    if(ReadLine != "#Init_description_meson_masses") crash("Error in input file. Expected: #Init_description_meson_masses, obtained: "+ReadLine);
    ReadInput>>ReadLine;
    if(ReadLine != "#Meson_to_Analyze:") crash("Error in input file. Expected: #Meson_to_Analyze, obtained: "+ReadLine);
    ReadInput >> Meson_to_analyze;
    ReadInput >> ReadLine;
    if(ReadLine != "#CURRENT_TYPE:") crash ("Error in input file. Expected: #CURRENT_TYPE, obtained: "+ReadLine);
    ReadInput >> CURRENT_TYPE;
    ReadInput >> ReadLine;
    if(ReadLine != "#Include_disconnected?:") crash("Error in input file. Expected: #Include_disconnected?, obtained: "+ReadLine);
    ReadInput >> IncludeDisconnected;
  }


  else if(Analysis_Mode =="Form_factors_Nissa") {
    //do nothing
  }

  else if(Analysis_Mode =="Form_factors_Nissa_3d") {
    //do nothing
  }

  else if(Analysis_Mode =="Bs_mumu_gamma") {
    //do nothing
  }

  else if(Analysis_Mode == "Semileptonic") {
    //do nothing
  }

  
    
  else if(Analysis_Mode=="Form_factors") {
    //do nothing  
  }

  else if(Analysis_Mode == "Axion_l7") {
    //do nothing
  }

  else if(Analysis_Mode == "Neutral_pi_TM") {
    //do nothing

  }

  else if(Analysis_Mode =="Eta_TM") {
    //do nothing

  }

  else if(Analysis_Mode=="VMD") {
    //do nothing
  }

  else if(Analysis_Mode=="GM2") {
    //do nothing

  }

  else if(Analysis_Mode=="R_ratio") {
    //do nothing

  }

  else if(Analysis_Mode=="Fake_spectral") {


    //do nothing

  }

  else if(Analysis_Mode=="tau_decay") {


    //do nothing
  

  }

  else if(Analysis_Mode=="tau_decay_strange") {


    //do nothing
  

  }

  else if(Analysis_Mode=="magnetic_susc") {

    //do nothing

  }

  else if(Analysis_Mode=="HVP") {

    //do nothing

  }

  else if(Analysis_Mode=="HVP_strange") {

    //do nothing

  }

  else if(Analysis_Mode=="RC_analysis") {

    //do nothing

  }

  else if(Analysis_Mode=="scale_setting") {

    //do nothing

  }

  else if(Analysis_Mode=="l7_Weinberg") {

    //do nothing

  }

  else if(Analysis_Mode=="weak_annihilation_inclusive") {

    //do nothing
  }

  
  else if(Analysis_Mode=="RC_WI_analysis") {

    //do nothing
  }

  else if(Analysis_Mode=="test") {

    //do nothing

  }

  
  
  else crash("Analysis_Mode: "+Analysis_Mode+" not found");
  
  ReadInput.close();
  Analysis_manager();
 
  return;

}


void MasterClass_analysis::Analysis_manager() {

  if(Analysis_Mode=="Form_factors") Compute_form_factors();

  if(Analysis_Mode=="Form_factors_Nissa") Compute_form_factors_Nissa();

  if(Analysis_Mode=="Form_factors_Nissa_3d") Get_radiative_form_factors_3d();

  if(Analysis_Mode=="Bs_mumu_gamma") Compute_Bs_mumu_gamma();

  if(Analysis_Mode=="Meson_masses") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode=="Meson_masses_twisted") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis_twisted(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode=="Meson_masses_ov_X") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis_ov_X(this->CURRENT_TYPE, this->IncludeDisconnected);
  }
  
  if(Analysis_Mode=="Meson_masses_twisted_adim") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis_twisted_adim(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode=="Meson_masses_twisted_ov_X") {
    if(Meson_to_analyze =="PI") Pion_mass_analysis_twisted_ov_X(this->CURRENT_TYPE, this->IncludeDisconnected);
  }

  if(Analysis_Mode == "Axion_l7") {
    Axion_l7_analysis();
  }

  if(Analysis_Mode == "Neutral_pi_TM") {
    Neutral_pi_TM();
  }

  if(Analysis_Mode == "Eta_TM") {
    Eta_TM();
  }

  if(Analysis_Mode == "VMD") {
    Regge();
  }

  if(Analysis_Mode == "GM2") {

    Gm2();

  }
  
  if(Analysis_Mode == "R_ratio") {

    R_ratio_analysis();
  }

  

  if(Analysis_Mode == "tau_decay") {

    tau_decay_analysis();
  }

  if(Analysis_Mode == "tau_decay_strange") {

    tau_decay_analysis_strange();
  }

  if(Analysis_Mode == "Semileptonic") {
    semileptonic_FF_analysis();
  }

  if(Analysis_Mode=="magnetic_susc") {
    Compute_magnetic_susc();
  }

  if(Analysis_Mode=="HVP") {
    HVP();
  }

  if(Analysis_Mode=="HVP_strange") {
    HVP_strange();
  }

  if(Analysis_Mode=="RC_analysis") {
    Perform_RC_analysis();
  }

  if(Analysis_Mode=="scale_setting") {
    Get_scale_setting();
  }

  if(Analysis_Mode=="l7_Weinberg") {
    l7_Weinberg();
  }

  if(Analysis_Mode=="weak_annihilation_inclusive") {

    compute_weak_annihilation();
  }

  if(Analysis_Mode=="RC_WI_analysis") {

    RC_WI_analysis();
  }

  
  if(Analysis_Mode=="test") {
    test();
  }
  
  return;


}
