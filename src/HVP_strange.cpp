#include "../include/HVP.h"


const double DTT = 0.5;
const double alpha = 1.0/137.035999;	     

using namespace std;








void HVP_strange() {





 
  

  int Njacks=50;
  bool UseJack=true;
  double qu= 2.0/3.0;
  double qd= -1.0/3.0;
  double fm_to_inv_Gev= 1.0/0.197327;

  //create directories
  boost::filesystem::create_directory("../data/HVP");
  boost::filesystem::create_directory("../data/HVP/Bounding");
  boost::filesystem::create_directory("../data/HVP/Corr");

  auto Sort_light_confs = [](string A, string B) {


    //return A<B;
			     
    int conf_length_A= A.length();
    int conf_length_B= B.length();

    int pos_a_slash=-1;
    int pos_b_slash=-1;
    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

    string A_bis= A.substr(pos_a_slash+1);
    string B_bis= B.substr(pos_b_slash+1);

   
			     
    string conf_num_A = A_bis.substr(0,4);
    string conf_num_B = B_bis.substr(0,4);
							       
		      
    string rA = A_bis.substr(A_bis.length()-5);
    string rB = B_bis.substr(B_bis.length()-5);
    if(rA.substr(0,1) == "r") { 
      int n1 = stoi(rA.substr(1,1));
      int n2 = stoi(rB.substr(1,1));
      if(rA == rB) {
	if(rA=="r0.h5" || rA=="r2.h5") return conf_num_A > conf_num_B;
	else if(rA=="r1.h5" || rA=="r3.h5") return conf_num_A < conf_num_B;
	else crash("stream not recognized");
      }
      else return n1<n2;
    }
    return A_bis<B_bis;
  };
  
  
 

  //############################################################################################
  //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
  GaussianMersenne GM(36551294);
  LatticeInfo a_info;
  distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
  distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack), ZV_E(UseJack);
  distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack), ZA_E(UseJack);
  double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
  double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err, ZV_E_ave, ZV_E_err;
  double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err, ZA_E_ave, ZA_E_err;
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp_FLAG;
  a_A_err= a_info.a_from_afp_FLAG_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp_FLAG;
  a_B_err= a_info.a_from_afp_FLAG_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp_FLAG;
  a_C_err= a_info.a_from_afp_FLAG_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp_FLAG;
  a_D_err= a_info.a_from_afp_FLAG_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp_FLAG;
  a_Z_err= a_info.a_from_afp_FLAG_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp_FLAG;
  a_E_err= a_info.a_from_afp_FLAG_err;
  ZA_E_ave = a_info.Za_WI_strange;
  ZA_E_err = a_info.Za_WI_strange_err;
  ZV_E_ave = a_info.Zv_WI_strange;
  ZV_E_err = a_info.Zv_WI_strange_err;

  for(int ijack=0;ijack<Njacks;ijack++) {

  a_A.distr.push_back( fm_to_inv_Gev*( a_A_ave + GM()*a_A_err*(1.0/sqrt(Njacks-1.0))));
  a_B.distr.push_back( fm_to_inv_Gev*( a_B_ave + GM()*a_B_err*(1.0/sqrt(Njacks-1.0))));
  a_C.distr.push_back( fm_to_inv_Gev*( a_C_ave + GM()*a_C_err*(1.0/sqrt(Njacks-1.0))));
  a_D.distr.push_back( fm_to_inv_Gev*( a_D_ave + GM()*a_D_err*(1.0/sqrt(Njacks-1.0))));
  a_Z.distr.push_back( fm_to_inv_Gev*( a_Z_ave + GM()*a_Z_err*(1.0/sqrt(Njacks-1.0))));
  a_E.distr.push_back( fm_to_inv_Gev*( a_E_ave + GM()*a_E_err*(1.0/sqrt(Njacks-1.0))));
  ZA_A.distr.push_back(  ZA_A_ave + GM()*ZA_A_err*(1.0/sqrt(Njacks -1.0)));
  ZV_A.distr.push_back(  ZV_A_ave + GM()*ZV_A_err*(1.0/sqrt(Njacks -1.0)));
  ZA_B.distr.push_back(  ZA_B_ave + GM()*ZA_B_err*(1.0/sqrt(Njacks -1.0)));
  ZV_B.distr.push_back(  ZV_B_ave + GM()*ZV_B_err*(1.0/sqrt(Njacks -1.0)));
  ZA_C.distr.push_back(  ZA_C_ave + GM()*ZA_C_err*(1.0/sqrt(Njacks -1.0)));
  ZV_C.distr.push_back(  ZV_C_ave + GM()*ZV_C_err*(1.0/sqrt(Njacks -1.0)));
  ZA_D.distr.push_back(  ZA_D_ave + GM()*ZA_D_err*(1.0/sqrt(Njacks -1.0)));
  ZV_D.distr.push_back(  ZV_D_ave + GM()*ZV_D_err*(1.0/sqrt(Njacks -1.0)));

  ZA_E.distr.push_back(  ZA_E_ave + GM()*ZA_E_err*(1.0/sqrt(Njacks -1.0)));
  ZV_E.distr.push_back(  ZV_E_ave + GM()*ZV_E_err*(1.0/sqrt(Njacks -1.0)));

  }
      

  //############################################################################################

















  bool Get_ASCII_strange= true;

    if(Get_ASCII_strange) {
    //read binary files
    boost::filesystem::create_directory("../HVP_strange");
    

    vector<string> Ens_T1({"B.72.64", "C.06.80", "C.06.112", "D.54.96", "E.44.112"});
    vector<string> Ens_TT1({"cB211b.072.64", "cC211a.06.80", "cC211a.06.112", "cD211a.054.96", "cE211a.044.112"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      vector<string> channels({"mix_s1_s1",  "mix_s2_s2"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../HVP_strange/"+channel);
	boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"TM_VKVK", "OS_VKVK"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../gm2_tau_rep_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
        size_t Nconfs, T, Nhits;
	bin_read(Nconfs, stream);
	bin_read(Nhits, stream);
	bin_read(T, stream);
	cout<<"channel: "<<channel<<endl;
	cout<<"Nconfs: "<<Nconfs<<endl;
	cout<<"T: "<<T<<" "<<T/2+1<<endl;
	cout<<"Nhits: "<<Nhits<<endl;
	for(size_t iconf=0;iconf<Nconfs;iconf++) {
	  vector<double> C(T/2+1);
	  for(size_t t=0;t<T/2+1;t++) bin_read(C[t], stream);
	  boost::filesystem::create_directory("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../HVP_strange/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
	  PrintCorr.precision(16);
	  PrintCorr<<"# "<<Corr_tags[id].substr(3,4)<<endl;
	  for(size_t t=0;t<(T/2+1);t++) PrintCorr<<C[t]<<endl;
	  if(Corr_tags[id].substr(3,4) == "VKTK" || Corr_tags[id].substr(3,4) == "TKVK") { for(size_t t=T/2+1; t<T;t++) PrintCorr<<-1*C[T-t]<<endl;   }
	  else  {for(size_t t=T/2+1; t<T;t++) PrintCorr<<C[T-t]<<endl;}
	  PrintCorr.close();
	
	}

	fclose(stream);

	}
	
      }
    }
    }

    

    auto Sort_easy = [](string A, string B) {

      int conf_length_A= A.length();
      int conf_length_B= B.length();
      
      int pos_a_slash=-1;
      int pos_b_slash=-1;
      for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
      for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;
      
      string A_bis= A.substr(pos_a_slash+1);
      string B_bis= B.substr(pos_b_slash+1);

      return atoi( A_bis.c_str()) < atoi( B_bis.c_str());
      
    };
   



    boost::filesystem::create_directory("../data/HVP_strange");

    data_t Vk_data_tm, Vk_data_OS;
    
    Vk_data_tm.Read("../HVP_strange/mix_s1_s1", "mes_contr_mix_s1_s1_TM_VKVK", "VKVK", Sort_easy);
    Vk_data_OS.Read("../HVP_strange/mix_s1_s1", "mes_contr_mix_s1_s1_OS_VKVK", "VKVK", Sort_easy);
    
    int Nens= Vk_data_tm.size;


    for(int iens=0;iens<Nens;iens++) {
      
      cout<<"Analyzing ensemble: "<<Vk_data_tm.Tag[iens]<<endl; 
      
      CorrAnalysis Corr(UseJack, Njacks,800);
      Corr.Nt = Vk_data_tm.nrows[iens];
     
      distr_t a_distr(UseJack);
      distr_t Zv(UseJack), Za(UseJack);

      if(Vk_data_tm.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; }
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
      else if(Vk_data_tm.Tag[iens].substr(1,1)=="E") {a_distr=a_E; Zv = ZV_E; Za = ZA_E;}
      else crash("Ensemble not found");
      
      LatticeInfo L_info;
      L_info.LatInfo_new_ens(Vk_data_tm.Tag[iens]);
      
      int L= L_info.L;

      
      
      distr_t_list Vk_tm_distr = 1e10*Za*Za*(pow(qd,2))*Corr.corr_t( Vk_data_tm.col(0)[iens] , "../data/HVP_strange/Corr/Vk_tm_"+Vk_data_tm.Tag[iens]+".dat");
      distr_t_list Vk_OS_distr = 1e10*Zv*Zv*(pow(qd,2))*Corr.corr_t( Vk_data_OS.col(0)[iens] , "../data/HVP_strange/Corr/Vk_OS_"+Vk_data_OS.Tag[iens]+".dat");
           
        
   
      distr_t amu_HVP_tm(UseJack,Njacks), amu_HVP_OS(UseJack,Njacks);

      auto K = [&](double Mv, double t, int size) -> double { return kernel_K(t, Mv);};

      distr_t_list Ker = distr_t_list::f_of_distr(K, a_distr , Corr.Nt/2);

      for(int t=1;t<Corr.Nt/2;t++) {

	amu_HVP_tm = amu_HVP_tm +  4.0*w(t,2)*pow(alpha,2)*Vk_tm_distr.distr_list[t]*Ker.distr_list[t];
	amu_HVP_OS = amu_HVP_OS +  4.0*w(t,2)*pow(alpha,2)*Vk_OS_distr.distr_list[t]*Ker.distr_list[t];

      }

    cout<<"#### "<<Vk_data_tm.Tag[iens]<<" ###"<<endl;
    cout<<"HVP tm: "<<amu_HVP_tm.ave()<<" +- "<<amu_HVP_tm.err()<<" stat. "<< (amu_HVP_tm.err()*100/amu_HVP_tm.ave())<<"%"<<endl;
    cout<<"HVP OS: "<<amu_HVP_OS.ave()<<" +- "<<amu_HVP_OS.err()<<" stat. "<< (amu_HVP_OS.err()*100/amu_HVP_OS.ave())<<"%"<<endl;
    cout<<"#######"<<endl;

         
    }


  return;

}
