#include "../include/weak_annihilation_inclusive.h"

using namespace std;


void compute_weak_annihilation() {

  int Njacks=59;
  bool UseJack=1;
  double fm_to_inv_Gev= 1.0/0.197327;

  boost::filesystem::create_directory("../data/WA_inclusive");

  auto Sort_confs = [](string A, string B) {

			   

			    int conf_length_A= A.length();
			    int conf_length_B= B.length();

			    int pos_a_slash=-1;
			    int pos_b_slash=-1;
			    for(int i=0;i<conf_length_A;i++) if(A.substr(i,1)=="/") pos_a_slash=i;
			    for(int j=0;j<conf_length_B;j++) if(B.substr(j,1)=="/") pos_b_slash=j;

			    string A_bis= A.substr(pos_a_slash+1);
			    string B_bis= B.substr(pos_b_slash+1);

			    //A_bis=A;
			    //B_bis=B;

			     
			    string conf_num_A = A_bis.substr(0,4);
			    string conf_num_B = B_bis.substr(0,4);
							       
		      
			    string rA = A_bis.substr(A_bis.length()-2);
			    string rB = B_bis.substr(B_bis.length()-2);
			    if(rA.substr(0,1) == "r") { 
			      int n1 = stoi(A_bis.substr(A_bis.length()-1));
			      int n2 = stoi(B_bis.substr(B_bis.length()-1));
			      if(rA == rB) {
			      if(rA=="r0" || rA=="r2") return conf_num_A > conf_num_B;
			      else if(rA=="r1" || rA=="r3") return conf_num_A < conf_num_B;
			      else crash("stream not recognized");
			      }
			      else return n1<n2;
			    }
			    return A_bis<B_bis;
  };
  
  

    vector<data_t> V_sou(4), A_sou(4);
    vector<data_t> V_sink(4), A_sink(4);

    data_t P5_sou, P5_sink;

    data_t P5_test;


    P5_test.Read("../test_weak_ann_mass", "mes_contr_H_B_H_H_S0_H", "P5P5", Sort_confs);

    int Nu= P5_test.Tag.size() ;

    for(int i=0;i<Nu;i++) {

      cout<<"Analyzing ensemble: "<<P5_test.Tag[i]<<endl; 

      CorrAnalysis Corr(UseJack, Njacks,800);
      Corr.Nt = P5_test.nrows[i];

      distr_t_list eff_mass = Corr.effective_mass_t( P5_test.col(0)[i], "../data/WA_inclusive/P5_test_"+P5_test.Tag[i] );


    }

    exit(-1);

    for(int i=0;i<4; i++) {
    
      V_sou[i].Read("../weak_annihilation_data", "mes_contr_H_B_PH_S0_H", "V"+to_string(i)+"P5", Sort_confs);
      A_sou[i].Read("../weak_annihilation_data", "mes_contr_H_B_PH_S0_H", "A"+to_string(i)+"P5", Sort_confs);

      V_sink[i].Read("../weak_annihilation_data", "mes_contr_H_B_PH_S0_H", "V"+to_string(i)+"P5", Sort_confs);
      A_sink[i].Read("../weak_annihilation_data", "mes_contr_H_B_PH_S0_H", "A"+to_string(i)+"P5", Sort_confs);

    }

    P5_sou.Read("../weak_annihilation_data", "mes_contr_H_B_H_H_S0_H" , "P5P5", Sort_confs);
    P5_sink.Read("../weak_annihilation_data", "mes_contr_H_B_H_H_S64_H", "P5P5", Sort_confs);


    int Nens=V_sou[0].Tag.size();



      //############################################################################################
    //generate fake jack_distr for lattice spacing a_A a_B, a_C, a_D and RENORMALIZATION CONSTANT
    GaussianMersenne GM(36551294);
    LatticeInfo a_info;
    distr_t a_A(UseJack), a_B(UseJack), a_C(UseJack), a_D(UseJack), a_Z(UseJack), a_E(UseJack);
    distr_t ZV_A(UseJack), ZV_B(UseJack), ZV_C(UseJack), ZV_D(UseJack);
    distr_t ZA_A(UseJack), ZA_B(UseJack), ZA_C(UseJack), ZA_D(UseJack);
    double a_A_ave, a_A_err, a_B_ave, a_B_err, a_C_ave, a_C_err, a_D_ave, a_D_err, a_Z_ave, a_Z_err, a_E_ave, a_E_err;
    double ZV_A_ave, ZV_A_err, ZV_B_ave, ZV_B_err, ZV_C_ave, ZV_C_err, ZV_D_ave, ZV_D_err;
    double ZA_A_ave, ZA_A_err, ZA_B_ave, ZA_B_err, ZA_C_ave, ZA_C_err, ZA_D_ave, ZA_D_err;
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
  
    }
      

  //############################################################################################






    boost::filesystem::create_directory("../data/WA_inclusive");

    

    for(int iens=0;iens<Nens;iens++) {

      cout<<"Analyzing ensemble: "<<P5_sou.Tag[iens]<<endl; 

      CorrAnalysis Corr(UseJack, Njacks,800);
      Corr.Nt = P5_sou.nrows[iens];
     
      distr_t a_distr(UseJack);
      distr_t Zv(UseJack), Za(UseJack);
      
      if(P5_sou.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; }
      else if(P5_sou.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C;}
      else if(P5_sou.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D;}
      else crash("Ensemble not found");

      LatticeInfo L_info;
      L_info.LatInfo_new_ens(P5_sou.Tag[iens]);
      
      //i
      VVfloat C_A1_sou = A_sink[1].col(0)[iens];
      VVfloat C_V1_sou = V_sink[1].col(0)[iens];
      VVfloat C_A1_sink = A_sink[1].col(0)[iens];
      VVfloat C_V1_sink = V_sink[1].col(0)[iens];
      VVfloat C_A2_sou = A_sink[2].col(0)[iens];
      VVfloat C_V2_sou = V_sink[2].col(0)[iens];
      VVfloat C_A2_sink = A_sink[2].col(0)[iens];
      VVfloat C_V2_sink = V_sink[2].col(0)[iens];

      //0
      VVfloat C_A0_sou = A_sink[0].col(0)[iens];
      VVfloat C_V0_sou = V_sink[0].col(0)[iens];
      VVfloat C_A0_sink = A_sink[0].col(0)[iens];
      VVfloat C_V0_sink = V_sink[0].col(0)[iens];
      //3
      VVfloat C_A3_sou = A_sink[3].col(0)[iens];
      VVfloat C_V3_sou = V_sink[3].col(0)[iens];
      VVfloat C_A3_sink = A_sink[3].col(0)[iens];
      VVfloat C_V3_sink = V_sink[3].col(0)[iens];
      

      //Y
      VVVfloat Y0_VV(Corr.Nt), Y1_VV(Corr.Nt), Y2_VV(Corr.Nt), Y3_VV(Corr.Nt);
      VVVfloat Y0_AA(Corr.Nt), Y1_AA(Corr.Nt), Y2_AA(Corr.Nt), Y3_AA(Corr.Nt);

      for(int i=0;i<Corr.Nt;i++) {
	Y0_VV[i].resize(Corr.Nt);
	Y1_VV[i].resize(Corr.Nt);
	Y2_VV[i].resize(Corr.Nt);
	Y3_VV[i].resize(Corr.Nt);

	Y0_AA[i].resize(Corr.Nt);
	Y1_AA[i].resize(Corr.Nt);
	Y2_AA[i].resize(Corr.Nt);
	Y3_AA[i].resize(Corr.Nt);
      }

      for(int i=0;i<Corr.Nt;i++) {
	for(int j=0;j<Corr.Nt;j++) {

	  Y0_VV[i][j] = summ_master( multiply_master( C_V1_sou[i], C_V1_sink[j] ), multiply_master( C_V2_sou[i], C_V2_sink[j])) ;
	  Y1_VV[i][j] = multiply_master( C_V0_sou[i], C_V0_sink[j] );
	  Y2_VV[i][j] = multiply_master( C_V3_sou[i], C_V3_sink[j] );
	  Y3_VV[i][j] = summ_master( multiply_master( C_V3_sou[i], C_V0_sink[j] ), multiply_master( C_V0_sou[i], C_V3_sink[j]));

	  Y0_AA[i][j] = summ_master( multiply_master( C_A1_sou[i], C_A1_sink[j] ), multiply_master( C_A2_sou[i], C_A2_sink[j])) ;
	  Y1_AA[i][j] = multiply_master( C_A0_sou[i], C_A0_sink[j] );
	  Y2_AA[i][j] = multiply_master( C_A3_sou[i], C_A3_sink[j] );
	  Y3_AA[i][j] = summ_master( multiply_master( C_A3_sou[i], C_A0_sink[j] ), multiply_master( C_A0_sou[i], C_A3_sink[j]));
	  

	}
      }
      
      
  
      distr_t_list P5_distr_sou(UseJack), P5_distr_sink(UseJack);

      P5_distr_sou = Corr.corr_t( P5_sou.col(0)[iens], "../data/WA_inclusive/P5_sou_"+P5_sou.Tag[iens] );
      P5_distr_sink = Corr.corr_t( P5_sink.col(0)[iens], "../data/WA_inclusive/P5_sink_"+P5_sou.Tag[iens] );



      

      
				  
      


    }


  



}
