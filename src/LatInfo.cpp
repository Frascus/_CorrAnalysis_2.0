#include "../include/LatInfo.h"


using namespace std;

const string Path_to_covariate_input_pars = "../Nf211_ETM_Ensemble_Info/CovMatrix_Branch";



void LatticeInfo::LatInfo(string S) {

  Za=0.0;
  Za_err=1.0;
  Za_M2=0.0;
  Za_M2_err= 1.0;
  double resc=1.0;


  if(S=="A") {
    Beta= 1.90;
    ZP= 0.529;
    ZP_err=0.007;
    ZS=0.747;
    ZS_err=0.012;
    a = 0.0885;
    a_err = 0.0036;
    ainv= 2.224;
    ainv_err= 0.068;
    Zm_fact=1.629;
    Zm_fact_err= 0.041;
    Za_fact=0.859;
    Za_fact_err=0.015;
    Zm_fact_M2= 1.637;
    Zm_fact_M2_err=0.014;
    Za_fact_M2= 0.990;
    Za_fact_M2_err=0.009;
    
    if(CURRENT_TYPE=="LOCAL") {
      Zv = 0.587;
      Zv_err=0.004;
      Zv_M2=0.608;
      Zv_M2_err=0.003;
      Za = 0.731*resc;
      Za_err = 0.008*resc;
      Za_M2 = 0.703*resc;
      Za_M2_err = 0.002*resc;
      
    }
    else if (CURRENT_TYPE=="CONSERVED") {Zv= 1.0; Zv_err=1.0; Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }
  else if(S=="B") {
    Beta= 1.95;
    ZP= 0.509;
    ZP_err=0.004;
    ZS=0.713;
    ZS_err=0.009;
    a = 0.0815;
    a_err = 0.0030;
    ainv= 2.416;
    ainv_err= 0.063;
    Zm_fact=1.514;
    Zm_fact_err= 0.033;
    Za_fact=0.873;
    Za_fact_err=0.013;
    Zm_fact_M2= 1.585;
    Zm_fact_M2_err=0.012;
    Za_fact_M2= 0.980;
    Za_fact_M2_err=0.008;
    if(CURRENT_TYPE=="LOCAL") {
      Zv=0.603;
      Zv_err=0.003;
      Zv_M2 = 0.614;
      Zv_M2_err = 0.002;
      Za = 0.737*resc;
      Za_err = 0.005*resc;
      Za_M2 = 0.714*resc;
      Za_M2_err = 0.002*resc;
      
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv= 1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
    
  }
  
  else if(S=="D") {
    Beta= 2.10;
    ZP= 0.516;
    ZP_err=0.02;
    ZS=0.700;
    ZS_err=0.006;
    a = 0.0619;
    a_err = 0.0018;
    ainv= 3.184;
    ainv_err= 0.059;
    Zm_fact=1.459;
    Zm_fact_err= 0.017;
    Za_fact=0.909;
    Za_fact_err=0.006;
    Zm_fact_M2= 1.462;
    Zm_fact_M2_err=0.006;
    Za_fact_M2= 0.958;
    Za_fact_M2_err=0.003;
    if(CURRENT_TYPE=="LOCAL") {
      Zv= 0.655 ;
      Zv_err=0.003;
      Zv_M2 = 0.657;
      Zv_M2_err = 0.002;
      Za = 0.762*resc;
      Za_err = 0.004*resc;
      Za_M2 = 0.752*resc;
      Za_M2_err = 0.002*resc;
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv=1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }

   else if(S=="a" || S=="b") {
    Beta= 2.10;
    ZP= 0.733207;
    ZP_err=0.00152823;
    //ZP=1.0;
    //ZP_err= 0.0000001;
    ZS=1.0;
    ZS_err=0.00000;
    a = 0.0619;
    a_err = 0.0018;
    if(CURRENT_TYPE=="LOCAL") {
      Zv= 0.655 ;
      Zv_err=0.003;
      Zv_M2 = 0.657;
      Zv_M2_err = 0.002;
      Za = 0.762;
      Za_err = 0.004;
      Za_M2 = 0.752;
      Za_M2_err = 0.002;
    }
    else if(CURRENT_TYPE=="CONSERVED") {Zv=1.0; Zv_err=1.0;  Zv_M2=1.0; Zv_M2_err=1.0;}
    else crash("CURRENT_TYPE: "+CURRENT_TYPE+" not yet implemented");
  }

  
  else crash("Simulation point: "+S+" is not defined");

  return;


}


void LatticeInfo::LatInfo_new_ens(string Tag) {

  if(Tag.substr(1,1)=="A") {
    a= 0.09471;          //real is 0.09471 let us change it by 1.5 sigma
    a_err= 0.00039;
    a_from_afp =  0.0907593; //0.0908026;
    a_from_afp_err = 0.00053917; //0.000535517;
    a_from_afp_FLAG= 0.0906883;
    a_from_afp_FLAG_err= 0.000510322;
    a_nucleon=0.09295;
    a_nucleon_err=0.00047;
    Za=0.724;
    Za_err=0.009;
    Za_WI = 0.695347;
    Za_WI_err= 0.00216267;
    Zv=0.6960;
    Zv_err=0.0007;
    Zv_WI =0.690716;
    Zv_WI_err = 8.78188e-05;
    Zv_WI_charm_extr= 0.601735;
    Zv_WI_charm_extr_err=0.000634;
    Za_WI_charm_extr = 0.66016;
    Za_WI_charm_extr_err = 0.0014;

    ZT_RI2 =  0.817432;  //0.8173;
    ZT_RI2_err = 2*0.00329; // sqrt(pow(0.0013,2) + pow(0.0008,2));

    ZT_RI2_M3= 0.8100;
    ZT_RI2_M3_err= sqrt( pow(0.0015,2) + pow(0.0008,2));
    

    
    if(Tag == "cA211a.53.24") {
      Za_WI_strange= 0.728393; // 0.728838 ;
      Za_WI_strange_err = 0.0018; // 0.00160579;
      Zv_WI_strange = 0.687002; // 0.68670;
      Zv_WI_strange_err = 0.00015; // 0.00014;
    }
    else if(Tag == "cA211a.40.24") {
      Za_WI_strange= 0.728393; // 0.728838 ;
      Za_WI_strange_err = 0.0018; // 0.00160579;
      Zv_WI_strange = 0.687002; // 0.68670;
      Zv_WI_strange_err = 0.00015; // 0.00014;
    }
    else if(Tag == "cA211ab.30.32") {
      Za_WI_strange= 0.728393; // 0.728838 ;
      Za_WI_strange_err = 0.0018; // 0.00160579;
      Zv_WI_strange = 0.687002; // 0.68670;
      Zv_WI_strange_err = 0.00015; // 0.00014;
    }

    else if(Tag == "cA211a.12.48") {
      Za_WI_strange= 0.728393; // 0.728838 ;
      Za_WI_strange_err = 0.0018; // 0.00160579;
      Zv_WI_strange = 0.687002; // 0.68670;
      Zv_WI_strange_err = 0.00015; // 0.00014;
    }

    Za_WI_FLAG= Za_WI_strange;
    Za_WI_FLAG_err= Za_WI_strange_err;

    Zv_WI_FLAG= Zv_WI_strange;
    Zv_WI_FLAG_err = Zv_WI_strange_err;
    
    ms_L=0.02050; ms_M=0.02300;
    ms_L_new=0.02050; ms_M_new=0.02300;
    mc_L= 0.26500; mc_M=0.29000; mc_H=0.30000;
  }

  //#######

  else if(Tag.substr(1,1)=="Z") {
    a= 0.08829;
    a_err= 0.00012;
    a_from_afp= 0.08829; //0.079616;
    a_from_afp_err = 0.00012;
    a_from_afp_FLAG= 0.08829; //0.079616;
    a_from_afp_FLAG_err = 0.00012;
    a_nucleon=0.08829;
    a_nucleon_err=0.00012;
    Za=0;
    Za_err=1;

    ZT_RI2 = 0;
    ZT_RI2_err = 1;

    ZT_RI2_M3= 0;
    ZT_RI2_M3_err= 1;
    
    if(Tag == "cZ211a.077.64") {
      Za_WI = 0 ;
      Za_WI_err = 1;
      Zv_WI = 0;
      Zv_WI_err = 1;

      //
      Za_WI_strange= 0;
      Za_WI_strange_err = 1;
      Zv_WI_strange = 0;
      Zv_WI_strange_err = 1;
    
    }

    Za_WI_FLAG= Za_WI_strange;
    Za_WI_FLAG_err= Za_WI_strange_err;

    Zv_WI_FLAG= Zv_WI_strange;
    Zv_WI_FLAG_err = Zv_WI_strange_err;
    
    Zv=0;
    Zv_err=1;
    ms_L=0.019; ms_M=0.021;
    ms_L_new=0.019; ms_M_new=0.021;
    mc_L= 0.21000; mc_M=0.23000; mc_H=0.25000;
  }

  else if(Tag.substr(1,1)=="B") {
    a= 0.08161;
    a_err= 0.00030;
    a_from_afp= 0.0795739; //0.079616;
    a_from_afp_err = 0.000132632; // 0.000127363;
    a_from_afp_FLAG= 0.0795131;
    a_from_afp_FLAG_err=3.52596e-05;
    a_nucleon=0.07975;
    a_nucleon_err=0.00032;
    Za=0.746;
    Za_err=0.008;

    ZT_RI2 = 0.835279; //0.8403;
    ZT_RI2_err = 2*0.00309962; //sqrt(pow(0.0017,2) + pow(0.0008,2));

    ZT_RI2_M3= 0.8261;
    ZT_RI2_M3_err= sqrt( pow(0.0023,2) + pow(0.0008,2));
    
    if(Tag == "cB211b.072.64") {
      Za_WI = 0.71403 ;
      Za_WI_err = 0.00077 ;
      Zv_WI = 0.709932 ;
      Zv_WI_err = 0.000007;

      //
      Za_WI_strange= 0.742844;
      Za_WI_strange_err = 0.00026;
      Zv_WI_strange = 0.706382;
      Zv_WI_strange_err = 2.4e-5;


      Za_WI_FLAG= 0.74300192;
      Za_WI_FLAG_err= 0.00020844752;

      Zv_WI_FLAG= 0.70637654;
      Zv_WI_FLAG_err= 1.9986831e-05;

      ms_L=0.019; ms_M=0.021;
      ms_L_new=0.018; ms_M_new=0.020;
      mc_L= 0.21000; mc_M=0.23000; mc_H=0.25000;
    
    }
    else if (Tag == "cB211b.072.96") {
      Za_WI = 0.71577 ;
      Za_WI_err = 0.00035 ;
      Zv_WI = 0.709950 ;
      Zv_WI_err = 0.000005;

      //
     
      Za_WI_strange= 0.742669;
      Za_WI_strange_err = 0.00015;
      Zv_WI_strange = 0.706406;
      Zv_WI_strange_err = 1.7e-5;

      Za_WI_FLAG= 0.74278317;
      Za_WI_FLAG_err= 0.00020476846;

      Zv_WI_FLAG= 0.70642655;
      Zv_WI_FLAG_err= 9.5013688e-06;
      
      ms_L=0.019; ms_M=0.021;
      ms_L_new=0.018; ms_M_new=0.019;
      mc_L= 0.21000; mc_M=0.23000; mc_H=0.25000;
   
    }
    Zv=0.7131;
    Zv_err=0.0006;
   
  }

  else if(Tag.substr(1,1)=="C") {
    a= 0.06942; 
    a_err= 0.00026;
    a_from_afp = 0.0682083;// 0.0682068;
    a_from_afp_err  = 0.000134938;// 0.000117345;
    a_from_afp_FLAG=0.0681569;
    a_from_afp_FLAG_err=8.2454e-05;
    a_nucleon=0.06860;
    a_nucleon_err=0.00020;
    Za=0.761;
    Za_err=0.008;
    Za_WI = 0.73803 ;
    Za_WI_err = 0.00047 ;
    Zv=0.7310;
    Zv_err=0.0005;
    Zv_WI = 0.728477 ;
    Zv_WI_err = 0.000005;

    //
    Za_WI_strange= 0.758355;
    Za_WI_strange_err = 0.00016;

    Zv_WI_strange = 0.725404;
    Zv_WI_strange_err = 1.9e-5;

    if(Tag== "cC211a.06.80") {

      Za_WI_FLAG= 0.75814062;
      Za_WI_FLAG_err= 0.00012782033;
      Zv_WI_FLAG= 0.72540536;
      Zv_WI_FLAG_err= 1.3640038e-05;

    }
    else { // cC211a.06.112

      Za_WI_FLAG= 0.75827928;
      Za_WI_FLAG_err= 0.00010801694;
      Zv_WI_FLAG= 0.72542129;
      Zv_WI_FLAG_err= 1.0065254e-05;

    }


    ZT_RI2 = 0.85619; // 0.8623;
    ZT_RI2_err = 2*0.0028443 ; // sqrt(pow(0.0019,2) + pow(0.0008,2));

    ZT_RI2_M3= 0.8545;
    ZT_RI2_M3_err = sqrt( pow(0.0019,2) + pow(0.0008,2));

   
    
    ms_L= 0.01600; ms_M= 0.01800;
    ms_L_new= 0.01600; ms_M_new= 0.01800;
    mc_L=0.17500; mc_M=0.19500; mc_H=0.21500;
  }

  else if(Tag.substr(1,1)=="D") {
    a= 0.0577;
    a_err=  0.0002;
    a_from_afp = 0.0569183;// 0.0569252;
    a_from_afp_err = 0.000115387;// 0.000103587;
    a_from_afp_FLAG=0.0568756;
    a_from_afp_FLAG_err=5.89089e-05;
    a_nucleon= 0.05777; //fake
    a_nucleon_err= 0.0002; //fake
    Za= 0.76217;
    Za_err=0.00024;
    Za_WI = 0.76134 ;
    Za_WI_err = 0.00027 ;
    Zv_WI = 0.746595 ;
    Zv_WI_err = 0.000005;
    Zv= 0.746599;
    Zv_err = 0.000005;

    //
    Za_WI_strange= 0.773944;
    Za_WI_strange_err = 0.00014;
    
    Zv_WI_strange = 0.744106;
    Zv_WI_strange_err = 1.2e-5;

    Za_WI_FLAG= 0.77366855;
    Za_WI_FLAG_err= 7.6206211e-05;

    Zv_WI_FLAG= 0.7441097;
    Zv_WI_FLAG_err= 6.9770149e-06;

    ZT_RI2 =  0.879164; //0.8872;
    ZT_RI2_err = 2*0.00246186; //sqrt(pow(0.0008,2) + pow(0.0008,2));

    ZT_RI2_M3= 0.8785;
    ZT_RI2_M3_err= sqrt( pow(0.0008,2) + pow(0.0008,2));

      
    ms_L= 0.014; ms_M= 0.015;
    ms_L_new= 0.013; ms_M_new= 0.014;
    mc_L=0.165; mc_M=0.175; mc_H=0.175; //mc_H is fake 



  }

  else if(Tag.substr(1,1)=="E") {
    a= 0;
    a_err=  0.0;
    a_from_afp = 0.0488471;
    a_from_afp_err = 5.34768e-05;
    a_from_afp_FLAG= 0.0489061;
    a_from_afp_FLAG_err = 5.99147e-05;
    a_nucleon= 0.0; //fake
    a_nucleon_err= 0.0; //fake
    Za= 0.0;
    Za_err=0.0;
    Za_WI = 0.0 ;
    Za_WI_err = 0.0 ;
    Zv_WI = 0.0 ;
    Zv_WI_err = 0.0;
    Zv= 0.0;
    Zv_err = 0.0;

    //
    Za_WI_strange= 0.7854224386568596;
    Za_WI_strange_err = 7.37626786405621e-05;
    
    Zv_WI_strange = 0.7582316794779995;
    Zv_WI_strange_err = 6.790392603635336e-06;


    Za_WI_FLAG= 0.78541808;
    Za_WI_FLAG_err= 7.244954e-05;
 
    Zv_WI_FLAG= 0.75823119;
    Zv_WI_FLAG_err= 5.0837341e-06; 
 

 

    ZT_RI2 = 0.0;
    ZT_RI2_err = 0.0;

    ZT_RI2_M3= 0.0;
    ZT_RI2_M3_err= 0.0;

      
    ms_L= 0.011; ms_M= 0.012;
    ms_L_new= 0.011; ms_M_new= 0.012;
    mc_L=0.165; mc_M=0.175; mc_H=0.175; //mc_H is fake 



  }

  

  else crash("In LatticeInfo::LatInfo_new_ens Ensemble: "+Tag+" not found while reading lattice spacing");

  if(Tag=="cA211a.53.24") {
    L=24; T=48; ml=0.00530; 
  }
  else if(Tag=="cA211a.40.24") {
    L=24; T=48; ml=0.00400; 
  }
  else if(Tag=="cA211ab.30.32") {
    L=32; T=64; ml=0.00300;
  }
  else if(Tag=="cA211a.12.48") {
    L=48; T=96; ml=0.00120;
  }
  else if(Tag=="cZ211a.077.64") {
    L=64; T=2*64; ml=0.00077;
  }
  else if(Tag=="cB211a.25.24") {
    L=24; T=48; ml=0.00250;
  }
  else if(Tag=="cB211a.25.32") {
    L=32; T=64; ml=0.00250;
  }
  else if(Tag=="cB211a.25.48") {
    L=48; T=96; ml=0.00250;
  }
  else if(Tag=="cB211a.14.64") {
    L=64; T=128; ml=0.00140; 
  }
  else if(Tag=="cB211b.072.64") {
    L=64; T=128; ml=0.00072;
  }
  else if(Tag=="cB211b.072.96") {
    L=96; T=192; ml=0.00072;
  }
  else if(Tag=="cC211a.20.48") {
    L=48; T=96; ml=0.00200;
  }
  else if(Tag=="cC211a.06.80") {
    L=80; T=160; ml=0.00060; 
  }
  else if(Tag=="cC211a.06.112") {
    L=112; T=2*112; ml=0.00060; 
  }
  
  else if(Tag=="cD211a.054.96") {
    L=96; T=192; ml=0.00054;
  }

  else if(Tag=="cE211a.044.112") {
    L=112; T=2*112; ml=0.00044;
  }
  else crash("In LatticeInfo::LatInfo_new_ens Ensemble: "+Tag+" not found");
  
   


  return;
}




void Read_pars_from_ensemble_tag(vector<string>& Ens_vec, Vfloat& m_lat, Vfloat& L, Vfloat& T) {


  for (auto &Ens: Ens_vec) {

    m_lat.push_back(0.0001*stod(Ens.substr(1, Ens.find(".") -1)));
    L.push_back((double)stoi(Ens.substr(Ens.find(".")+1,Ens.find("_") -1 - Ens.find("."))));
    T.push_back((double)stoi(Ens.substr(Ens.find("_")+1, string::npos)));
  }

  return;
}
void Read_pars_from_ensemble_tag(string Ens, double &m_lat, double &L, double &T) {


   m_lat =0.0001*stod(Ens.substr(1, Ens.find(".") -1));
   L =  (double)stoi(Ens.substr(Ens.find(".")+1,Ens.find("_") -1 - Ens.find(".")));
   T = (double)stoi(Ens.substr(Ens.find("_")+1, string::npos));


   return;


}


void ReadBranch(int k, Eigen::MatrixXd& CovMatrixInput, Eigen::VectorXd& Ave_input_parameters ) {

  //read averages from file
  CovMatrixInput.resize(9,9);
  Ave_input_parameters.resize(9);
  ifstream ReadFromFile(Path_to_covariate_input_pars+to_string(k)+".txt");
  double m_l, B0, f0,ainv_temp1, ainv_temp2, ainv_temp3, Zp_temp1, Zp_temp2, Zp_temp3;
  ReadFromFile>>m_l>>B0>>f0>>ainv_temp1>>ainv_temp2>>ainv_temp3>>Zp_temp1>>Zp_temp2>>Zp_temp3;
  Ave_input_parameters<<m_l,B0,f0,ainv_temp1,ainv_temp2,ainv_temp3,Zp_temp1,Zp_temp2,Zp_temp3;
 


  string CommentLine="";
  while(CommentLine=="") getline(ReadFromFile, CommentLine);
  string Expected= "# Cov Matrix";
  if(CommentLine != Expected) crash("Error while reading CovMatrix, branch nr. "+to_string(k)+"\n Expected: '# Cov Matrix' \t Obtained: "+CommentLine);
  //read Covariance Matrix
  for(int i=0; i<9;i++) {
    for(int j=0;j<9;j++) {
      double temp;
      ReadFromFile >> temp;
      CovMatrixInput(i,j) = temp;
    }
  }

  double stop;
  ReadFromFile >> stop;
  if(!ReadFromFile.eof()) crash("Error: after reading branch nr. "+to_string(k)+" File is !eof");


  ReadFromFile.close();
  cout<<"LOADED COVARIANCE MATRIX FOR BRANCH nr. "<<k<<endl;
  cout<<Ave_input_parameters<<endl;
  cout<<CovMatrixInput<<endl;
 
  
  return;
}

