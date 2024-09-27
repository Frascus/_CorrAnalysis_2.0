#include "../include/l7_Weinberg.h"


using namespace std;





void l7_Weinberg() {



  int Njacks=50;
  bool UseJack=true;
  double fm_to_inv_Gev= 1.0/0.197327;

  
bool Get_ASCII= true;

    if(Get_ASCII) {
    //read binary files
    boost::filesystem::create_directory("../l7_Weinberg");
    

    vector<string> Ens_T1({"B.072.64", "C.06.80", "D.054.96"});
    vector<string> Ens_TT1({"cB211b.072.64", "cC211a.06.80", "cD211a.054.96"});

    for( int it=0; it<(signed)Ens_T1.size(); it++) {

      cout<<"Ens: "<<Ens_T1[it]<<endl;

      vector<string> channels({"ll"});

      for(auto &channel : channels) {
	boost::filesystem::create_directory("../l7_Weinberg/"+channel);
	boost::filesystem::create_directory("../l7_Weinberg/"+channel+"/"+Ens_TT1[it]);
      }
      //read binary
      vector<string> Corr_tags({"OS_P5P5", "TM_S0S0"});

          
      for(int id=0; id<(signed)Corr_tags.size(); id++) {

	cout<<"Corr: "<<Corr_tags[id]<<endl;
	for( auto &channel: channels) {

	FILE *stream = fopen( ("../l7_data_bin/"+Ens_T1[it]+"/"+channel+"_"+Corr_tags[id]).c_str(), "rb");
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
	  boost::filesystem::create_directory("../l7_Weinberg/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf));
	  ofstream PrintCorr("../l7_Weinberg/"+channel+"/"+Ens_TT1[it]+"/"+to_string(iconf)+"/mes_contr_"+channel+"_"+Corr_tags[id]);
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

			    //A_bis=A;
			    //B_bis=B;

			     
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


   


    data_t P5P5, S0S0, disc_P5P5;


    P5P5.Read("../l7_Weinberg/ll", "mes_contr_ll_OS_P5P5", "P5P5", Sort_easy);
    S0S0.Read("../l7_Weinberg/ll", "mes_contr_ll_TM_S0S0", "S0S0", Sort_easy);
    disc_P5P5.Read("../magnetic_susc_disco/disco_P5", "P5P5.txt", "P5P5", Sort_light_confs);





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
  
  distr_t ZpZs_A(UseJack), ZpZs_B(UseJack), ZpZs_C(UseJack), ZpZs_D(UseJack);
  a_info.LatInfo_new_ens("cA211a.53.24");
  a_A_ave= a_info.a_from_afp;
  a_A_err= a_info.a_from_afp_err;
  ZA_A_ave = a_info.Za_WI_strange;
  ZA_A_err = a_info.Za_WI_strange_err;
  ZV_A_ave = a_info.Zv_WI_strange;
  ZV_A_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cB211b.072.64");
  a_B_ave= a_info.a_from_afp;
  a_B_err= a_info.a_from_afp_err;
  ZA_B_ave = a_info.Za_WI_strange;
  ZA_B_err = a_info.Za_WI_strange_err;
  ZV_B_ave = a_info.Zv_WI_strange;
  ZV_B_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cC211a.06.80");
  a_C_ave= a_info.a_from_afp;
  a_C_err= a_info.a_from_afp_err;
  ZA_C_ave = a_info.Za_WI_strange;
  ZA_C_err = a_info.Za_WI_strange_err;
  ZV_C_ave = a_info.Zv_WI_strange;
  ZV_C_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cD211a.054.96");
  a_D_ave= a_info.a_from_afp;
  a_D_err= a_info.a_from_afp_err;
  ZA_D_ave = a_info.Za_WI_strange;
  ZA_D_err = a_info.Za_WI_strange_err;
  ZV_D_ave = a_info.Zv_WI_strange;
  ZV_D_err = a_info.Zv_WI_strange_err;
  a_info.LatInfo_new_ens("cZ211a.077.64");
  a_Z_ave= a_info.a_from_afp;
  a_Z_err= a_info.a_from_afp_err;
  a_info.LatInfo_new_ens("cE211a.044.112");
  a_E_ave= a_info.a_from_afp;
  a_E_err= a_info.a_from_afp_err;

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

  ZpZs_A.distr.push_back( 0.7517 + GM()*( 0.0029  )*(1.0/sqrt(Njacks-1.0)));
  ZpZs_B.distr.push_back( 0.79066 + GM()*( 0.00023   )*(1.0/sqrt(Njacks-1.0)));
  ZpZs_C.distr.push_back( 0.82308 + GM()*( 0.00023  )*(1.0/sqrt(Njacks-1.0)));
  ZpZs_D.distr.push_back( 0.85095 + GM()*( 0.00018  )*(1.0/sqrt(Njacks-1.0)));
  
 

  

  }
      

  //############################################################################################


  int Nens= P5P5.size;


  for(int iens=0;iens<Nens;iens++) {

    cout<<"Analyzing ensemble: "<<P5P5.Tag[iens]<<endl; 

    CorrAnalysis Corr(UseJack, Njacks,800);
    Corr.Nt = P5P5.nrows[iens];
     
    distr_t a_distr(UseJack);
    distr_t Zv(UseJack), Za(UseJack);
    distr_t ZpZs(UseJack);

    if(P5P5.Tag[iens].substr(1,1)=="B") {a_distr=a_B; Zv = ZV_B; Za = ZA_B; ZpZs= ZpZs_B;}
    else if(P5P5.Tag[iens].substr(1,1)=="C") {a_distr=a_C; Zv = ZV_C; Za = ZA_C; ZpZs= ZpZs_C;}
    else if(P5P5.Tag[iens].substr(1,1)=="D") {a_distr=a_D; Zv = ZV_D; Za = ZA_D; ZpZs=ZpZs_D;}
    else crash("Ensemble not found");

    LatticeInfo L_info;
    L_info.LatInfo_new_ens(P5P5.Tag[iens]);

    int L= L_info.L;

   
    boost::filesystem::create_directory("../data/l7_Weinberg");

    distr_t_list P5P5_distr = Corr.corr_t( P5P5.col(0)[iens], "../data/l7_Weinberg/P5P5_conn_"+P5P5.Tag[iens]+".dat");
    distr_t_list P5P5_disc_distr = Corr.corr_t( disc_P5P5.col(0)[iens], "../data/l7_Weinberg/P5P5_disc_"+P5P5.Tag[iens]+".dat");
    distr_t_list S0S0_distr = Corr.corr_t( S0S0.col(0)[iens], "../data/l7_Weinberg/S0S0_"+P5P5.Tag[iens]+".dat");

    distr_t_list C = S0S0_distr - P5P5_distr +  2*P5P5_disc_distr;


    distr_t_list P5_full = P5P5_distr -2*P5P5_disc_distr;

    Print_To_File({}, {P5_full.ave(), P5_full.err()}, "../data/l7_Weinberg/P5P5_"+P5P5.Tag[iens]+".dat.t", "", "");

    Print_To_File({}, {C.ave(), C.err()}, "../data/l7_Weinberg/diff_c_"+P5P5.Tag[iens]+".dat.t", "", "");


    if(P5P5.Tag[iens].substr(1,12)=="B211b.072.96") {Corr.Tmin=30; Corr.Tmax=70;}
    else if(P5P5.Tag[iens].substr(1,12)=="B211b.072.64") { Corr.Tmin=27; Corr.Tmax=50;}
    else if(P5P5.Tag[iens].substr(1,1)=="C") {Corr.Tmin=40; Corr.Tmax=60;}
    else if(P5P5.Tag[iens].substr(1,1)=="D") {Corr.Tmin=41; Corr.Tmax=80;}


    distr_t Mpi = Corr.Fit_distr( Corr.effective_mass_t(P5P5_distr, ""));
    double ml= L_info.ml;

    distr_t F= -2*(ml*ml/POW_D(Mpi,4))*POW_D(1.0/ZpZs,2);

    distr_t_list l7_p_sum(UseJack);

    l7_p_sum.distr_list.push_back( 0.0*Get_id_jack_distr(Njacks));

    int t0=1;
    
    for(int t=t0;t<Corr.Nt/2;t++) { l7_p_sum.distr_list.push_back( l7_p_sum.distr_list[t-t0] + w(t,3)*F*C[t]);}

    distr_t l7= l7_p_sum[l7_p_sum.size()-1];

    cout<<"l7["<<P5P5.Tag[iens]<<"]: "<<l7.ave()<<" +- "<<l7.err()<<endl;
    cout<<"F: "<<F.ave()<<endl;

    Print_To_File({},{l7_p_sum.ave(), l7_p_sum.err()}, "../data/l7_Weinberg/l7_psum_"+P5P5.Tag[iens]+".dat", "", "");
    
    
      

      

   
       
  }

  
    
  return;



}
