#include "../include/Meson_mass_extrapolation.h"

const bool verb=true;
using namespace std;


class Ov {

public:
  Ov()  {}
  
  double Y,Y_err, X, X_phys;
};

class fpar {

public:
  fpar() {}
  fpar(const Vfloat &par) {
    if((signed)par.size() != 3) crash("In Obs_extrapolation_meson_mass, class fpar has size different from three.");
    a0=par[0];
    D=par[1];
    D2= par[2];
  }

  double a0,D, D2;

};



distr_t Obs_extrapolation_meson_mass( vector<distr_t> &Meas, vector<distr_t> &Op, double Op_phys, string Print_dir, string Tag, bool UseJack, string EXTRAPOLATION_MODE) {

  //Get number of jackknives
  double Njacks= Meas[0].size();
  double Nmeas = Meas.size();
  if(Nmeas != Op.size()) crash("Meas and Op size do not coincide");
  cout<<"Njacks: "<<Njacks<<endl;
  cout<<"Nmeas: "<<Nmeas<<endl;
  cout<<"Print directory: "<<Print_dir<<endl;
  cout<<"Tag: "<<Tag<<endl;
  cout<<"##########"<<endl;
  cout<<"Meson obs values (Exp: "<<Op_phys<<" )"<<endl;
  for(int imeas=0;imeas<Nmeas;imeas++) cout<<Op[imeas].ave()<<" +- "<<Op[imeas].err()<<" ; "<<Meas[imeas].ave()<<" +- "<<Meas[imeas].err()<<endl;
  cout<<"##########"<<endl;
  distr_t a0_fit(UseJack), D_fit(UseJack), D2_fit(UseJack);

  if(EXTRAPOLATION_MODE =="FIT") {
  bootstrap_fit<fpar,Ov> bf(Njacks);
  bf.set_warmup_lev(2); //sets warmup

  bf.Set_number_of_measurements(Meas.size());
  bf.Set_verbosity(verb);
  double a0_guess =0;
  for(auto &ms: Meas) a0_guess += ms.ave()/Nmeas;
  bf.Add_par("a0", a0_guess, 1e-2*a0_guess);
  bf.Add_par("D", -1.0e-1, 1e-3);
  bf.Add_par("D2", -1.0e-1, 1e-3);
  if(Nmeas <= 3 || Nmeas >= 7) bf.Fix_par("D2", 0.0);

  bf.ansatz =  [&](const fpar &p, const Ov &ip) -> double {
		 return p.a0*(1.0) + (p.D*(ip.X - Op_phys)) +  p.D2*pow((ip.X - Op_phys),2);
	       };

  bf.measurement = [&](const fpar &p, const Ov &ip) ->double {
		     return ip.Y;
		   };
  bf.error = [&](const fpar &p, const Ov &ip) -> double {
	       return ip.Y_err;
	     };


  //meas and errs

  
  //fill the data
  vector<vector<Ov>> data(Njacks);
  for(auto &data_ij: data) data_ij.resize(Nmeas);
  //allocate space for output result
  boot_fit_data<fpar> Bt_fit;

  for(auto &data_iboot: data) data_iboot.resize(Nmeas);
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int i=0;i<Nmeas;i++) {
      data[ijack][i].Y = Meas[i].distr[ijack];
      data[ijack][i].Y_err = Meas[i].err();
      data[ijack][i].X = Op[i].distr[ijack];
    } 
  }


    
  //append
  bf.Append_to_input_par(data);

  //fit
  Bt_fit= bf.Perform_bootstrap_fit();
  Bt_fit.ch2_ave();


  //retrieve fit parameter a0

  

  for(int ijack=0;ijack<Njacks;ijack++) {a0_fit.distr.push_back( Bt_fit.par[ijack].a0); D_fit.distr.push_back(Bt_fit.par[ijack].D); D2_fit.distr.push_back( Bt_fit.par[ijack].D2);}


  if(Print_dir != "") {
    boost::filesystem::create_directory(Print_dir+"/extr");
    vector<double> X, Xerr, Y,Yerr, Extr_val;
    for(int imeas=0;imeas<Nmeas;imeas++) {
      X.push_back( Op[imeas].ave());
      Xerr.push_back( Op[imeas].err());
      Y.push_back( Meas[imeas].ave());
      Yerr.push_back( Meas[imeas].err());
      Extr_val.push_back(0);
    }

    double der_val = (0.01*Op_phys*D_fit/a0_fit).ave();
    double der_err = (0.01*Op_phys*D_fit/a0_fit).err();
    X.push_back(Op_phys);
    Xerr.push_back(0.0);
    Y.push_back(a0_fit.ave());
    Yerr.push_back(a0_fit.err());
    Extr_val.push_back(1);

    Print_To_File({}, {X,Xerr, Y, Yerr, Extr_val}, Print_dir+"/extr/"+Tag+"_fit_extr.meas", "", "# X Xerr Y Yerr Extr_id,  Derivative: "+to_string_with_precision(der_val, 6)+" +- "+to_string_with_precision(der_err,6));

    //print fitting function
    Vfloat ave_ops;
    for(auto & op: Op) ave_ops.push_back( op.ave());
    int Npoints = 300;
    double xmin= min(0.5*Op_phys, 0.5*(*min_element(ave_ops.begin(), ave_ops.end())));
    double xmax= max(1.5*Op_phys, 1.5*(*max_element(ave_ops.begin(), ave_ops.end())));
    Vfloat Xpoints;
    distr_t_list Fit_func;
    for(int ipoint=0;ipoint<Npoints;ipoint++) {
      double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
      Fit_func.distr_list.push_back(a0_fit*(1.0) + (D_fit*(xp-Op_phys)) + D2_fit*pow( xp- Op_phys,2));
      Xpoints.push_back(xp);
    }

    Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},Print_dir+"/extr/"+Tag+"_fit_extr.fit_func", "", "");
    


  }


  }
  else if (EXTRAPOLATION_MODE == "SPLINE" ) {

    
    for(int ijack=0;ijack<Njacks;ijack++) {
      //find the quadratic function that passes through the points if Npoints = 3, or the linear function that passes through the points if Npoints = 2;
      if(Nmeas==2) {
	double y1 = Meas[0].distr[ijack];
	double y2 = Meas[1].distr[ijack];
	double Dx1 = Op[0].distr[ijack] - Op_phys;
	double Dx2 = Op[1].distr[ijack] -Op_phys;
	a0_fit.distr.push_back(  (y1*Dx2 - y2*Dx1)/(Dx2-Dx1));
	D_fit.distr.push_back( -y1/(Dx2-Dx1) +y2/(Dx2-Dx1) );
      }
      else if(Nmeas==3) {
	double y1 = Meas[0].distr[ijack];
	double y2 = Meas[1].distr[ijack];
	double y3 = Meas[2].distr[ijack];
	double Dx1 = Op[0].distr[ijack] - Op_phys;
	double Dx2 = Op[1].distr[ijack] -Op_phys;
	double Dx3 = Op[2].distr[ijack] -Op_phys;

	double n1 = y1*Dx2*Dx3/( (Dx1-Dx2)*(Dx1-Dx3));
	double n2 = -y2*Dx1*Dx3/( (Dx1-Dx2)*(Dx2-Dx3));
	double n3 = y3*Dx1*Dx2/( (Dx1-Dx3)*(Dx2-Dx3));

	double prod= (Dx1-Dx2)*(Dx1-Dx3)*(Dx2-Dx3);

	a0_fit.distr.push_back( n1+n2+n3);
	D_fit.distr.push_back( (pow(Dx3,2)*(y1-y2) + pow(Dx1,2)*(y2-y3) + pow(Dx2,2)*(y3-y1))/prod     );
	D2_fit.distr.push_back( (Dx3*(y2-y1) + Dx2*(y1-y3) + Dx1*(y3-y2))/prod   );




      }
      else crash("In Meson_mass_extrapolation.cpp: SPLINATOR called with Nmeas: "+to_string(Nmeas));
    }

    //print result

    
    if(Print_dir != "") {
      boost::filesystem::create_directory(Print_dir+"/extr");
      vector<double> X, Xerr, Y,Yerr, Extr_val;
      for(int imeas=0;imeas<Nmeas;imeas++) {
	X.push_back( Op[imeas].ave());
	Xerr.push_back( Op[imeas].err());
	Y.push_back( Meas[imeas].ave());
	Yerr.push_back( Meas[imeas].err());
	Extr_val.push_back(0);
      }
      
      X.push_back(Op_phys);
      Xerr.push_back(0.0);
      Y.push_back(a0_fit.ave());
      Yerr.push_back(a0_fit.err());
      Extr_val.push_back(1);

      double der_val= (0.01*Op_phys*D_fit/a0_fit).ave();
      double der_err= (0.01*Op_phys*D_fit/a0_fit).err();
      
      Print_To_File({}, {X,Xerr, Y, Yerr, Extr_val}, Print_dir+"/extr/"+Tag+"_spline_extr.meas", "", "# X Xerr Y Yerr Extr_id, derivative: "+to_string_with_precision(der_val,6)+" +- "+to_string_with_precision(der_err,6));
      
      //print fitting function
      Vfloat ave_ops;
      for(auto & op: Op) ave_ops.push_back( op.ave());
      int Npoints = 300;
      double xmin= min(0.5*Op_phys, 0.5*(*min_element(ave_ops.begin(), ave_ops.end())));
      double xmax= max(1.5*Op_phys, 1.5*(*max_element(ave_ops.begin(), ave_ops.end())));
      Vfloat Xpoints;
      distr_t_list Fit_func;
      for(int ipoint=0;ipoint<Npoints;ipoint++) {
	double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
	if(Nmeas==2) { Fit_func.distr_list.push_back(a0_fit + D_fit*(xp-Op_phys)) ;}
	else if(Nmeas==3)  { Fit_func.distr_list.push_back(a0_fit + D_fit*(xp-Op_phys) + D2_fit*(xp-Op_phys)*(xp-Op_phys)  );}
	else crash("In Meson_mass_extrapolation.cpp: SPLINATOR called with Nmeas: "+to_string(Nmeas));
	Xpoints.push_back(xp);
      }

    Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},Print_dir+"/extr/"+Tag+"_spline_extr.fit_func", "", "");
    


    }

    
  }

  return a0_fit;
  
}




distr_t Obs_extrapolation_meson_mass( vector<distr_t> &Meas, vector<distr_t> &Op, const distr_t  &Op_phys, string Print_dir, string Tag, bool UseJack, string EXTRAPOLATION_MODE) {

  //Get number of jackknives
  double Njacks= Meas[0].size();
  double Nmeas = Meas.size();
  if(Nmeas != Op.size()) crash("Meas and Op size do not coincide");
  cout<<"Njacks: "<<Njacks<<endl;
  cout<<"Nmeas: "<<Nmeas<<endl;
  cout<<"Print directory: "<<Print_dir<<endl;
  cout<<"Tag: "<<Tag<<endl;
  cout<<"##########"<<endl;
  cout<<"Meson obs values (Exp: "<<Op_phys.ave()<<" )"<<endl;
  for(int imeas=0;imeas<Nmeas;imeas++) cout<<Op[imeas].ave()<<" +- "<<Op[imeas].err()<<" ; "<<Meas[imeas].ave()<<" +- "<<Meas[imeas].err()<<endl;
  cout<<"##########"<<endl;
  distr_t a0_fit(UseJack), D_fit(UseJack), D2_fit(UseJack);

  if(EXTRAPOLATION_MODE =="FIT") {
  bootstrap_fit<fpar,Ov> bf(Njacks);
  bf.set_warmup_lev(2); //sets warmup

  bf.Set_number_of_measurements(Meas.size());
  bf.Set_verbosity(verb);
  double a0_guess =0;
  for(auto &ms: Meas) a0_guess += ms.ave()/Nmeas;
  bf.Add_par("a0", a0_guess, 1e-2*a0_guess);
  bf.Add_par("D", -1.0e-1, 1e-3);
  bf.Add_par("D2", -1.0e-1, 1e-3);
  if(Nmeas <= 3 || Nmeas > 7) bf.Fix_par("D2", 0.0);

  bf.ansatz =  [&](const fpar &p, const Ov &ip) -> double {
		 return p.a0*(1.0) + p.D*(ip.X - ip.X_phys) + p.D2*pow( ip.X - ip.X_phys, 2);
	       };

  bf.measurement = [&](const fpar &p, const Ov &ip) ->double {
		     return ip.Y;
		   };
  bf.error = [&](const fpar &p, const Ov &ip) -> double {
	       return ip.Y_err;
	     };


  //meas and errs

  
  //fill the data
  vector<vector<Ov>> data(Njacks);
  for(auto &data_ij: data) data_ij.resize(Nmeas);
  //allocate space for output result
  boot_fit_data<fpar> Bt_fit;

  for(auto &data_iboot: data) data_iboot.resize(Nmeas);
  
  for(int ijack=0;ijack<Njacks;ijack++) {
    for(int i=0;i<Nmeas;i++) {
      data[ijack][i].Y = Meas[i].distr[ijack];
      data[ijack][i].Y_err = Meas[i].err();
      data[ijack][i].X = Op[i].distr[ijack];
      data[ijack][i].X_phys = Op_phys.distr[ijack];
    } 
  }


    
  //append
  bf.Append_to_input_par(data);

  //fit
  Bt_fit= bf.Perform_bootstrap_fit();
  Bt_fit.ch2_ave();


  //retrieve fit parameter a0

  

  for(int ijack=0;ijack<Njacks;ijack++) {a0_fit.distr.push_back( Bt_fit.par[ijack].a0); D_fit.distr.push_back(Bt_fit.par[ijack].D); D2_fit.distr.push_back( Bt_fit.par[ijack].D2);}


  if(Print_dir != "") {
    boost::filesystem::create_directory(Print_dir+"/extr");
    vector<double> X, Xerr, Y,Yerr, Extr_val;
    for(int imeas=0;imeas<Nmeas;imeas++) {
      X.push_back( Op[imeas].ave());
      Xerr.push_back( Op[imeas].err());
      Y.push_back( Meas[imeas].ave());
      Yerr.push_back( Meas[imeas].err());
      Extr_val.push_back(0);
    }
   
    X.push_back(Op_phys.ave());
    Xerr.push_back(Op_phys.err());
    Y.push_back(a0_fit.ave());
    Yerr.push_back(a0_fit.err());
    Extr_val.push_back(1);

    double der_val =(0.01*Op_phys*D_fit/a0_fit).ave();
    double der_err= (0.01*Op_phys*D_fit/a0_fit).err();

    Print_To_File({}, {X,Xerr, Y, Yerr, Extr_val}, Print_dir+"/extr/"+Tag+"_fit_extr.meas", "", "# X Xerr Y Yerr Extr_id, derivative: "+to_string_with_precision(der_val,6)+" +- "+to_string_with_precision(der_err,6));

    //print fitting function
    Vfloat ave_ops;
    for(auto & op: Op) ave_ops.push_back( op.ave());
    int Npoints = 300;
    double xmin= min(0.5*Op_phys.ave(), 0.5*(*min_element(ave_ops.begin(), ave_ops.end())));
    double xmax= max(1.5*Op_phys.ave(), 1.5*(*max_element(ave_ops.begin(), ave_ops.end())));
    Vfloat Xpoints;
    distr_t_list Fit_func;
    for(int ipoint=0;ipoint<Npoints;ipoint++) {
      double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
      Fit_func.distr_list.push_back(a0_fit*(1.0) + D_fit*(xp-Op_phys) + D2_fit*(xp - Op_phys)*(xp - Op_phys));
      Xpoints.push_back(xp);
    }

    Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},Print_dir+"/extr/"+Tag+"_fit_extr.fit_func", "", "");
    


  }


  }
  else if (EXTRAPOLATION_MODE == "SPLINE" ) {

    
    for(int ijack=0;ijack<Njacks;ijack++) {
      //find the quadratic function that passes through the points if Npoints = 3, or the linear function that passes through the points if Npoints = 2;
      if(Nmeas==2) {
	double y1 = Meas[0].distr[ijack];
	double y2 = Meas[1].distr[ijack];
	double Dx1 = Op[0].distr[ijack] - Op_phys.distr[ijack];
	double Dx2 = Op[1].distr[ijack] -Op_phys.distr[ijack];
	a0_fit.distr.push_back(  (y1*Dx2 - y2*Dx1)/(Dx2-Dx1));
	D_fit.distr.push_back( -y1/(Dx2-Dx1) +y2/(Dx2-Dx1) );
      }
      else if(Nmeas==3) {
	double y1 = Meas[0].distr[ijack];
	double y2 = Meas[1].distr[ijack];
	double y3 = Meas[2].distr[ijack];
	double Dx1 = Op[0].distr[ijack] - Op_phys.distr[ijack];
	double Dx2 = Op[1].distr[ijack] -Op_phys.distr[ijack];
	double Dx3 = Op[2].distr[ijack] -Op_phys.distr[ijack];

	double n1 = y1*Dx2*Dx3/( (Dx1-Dx2)*(Dx1-Dx3));
	double n2 = -y2*Dx1*Dx3/( (Dx1-Dx2)*(Dx2-Dx3));
	double n3 = y3*Dx1*Dx2/( (Dx1-Dx3)*(Dx2-Dx3));

	double prod= (Dx1-Dx2)*(Dx1-Dx3)*(Dx2-Dx3);

	a0_fit.distr.push_back( n1+n2+n3);
	D_fit.distr.push_back( (pow(Dx3,2)*(y1-y2) + pow(Dx1,2)*(y2-y3) + pow(Dx2,2)*(y3-y1))/prod     );
	D2_fit.distr.push_back( (Dx3*(y2-y1) + Dx2*(y1-y3) + Dx1*(y3-y2))/prod   );


      }
      else crash("In Meson_mass_extrapolation.cpp: SPLINATOR called with Nmeas: "+to_string(Nmeas));
    }

    //print result

    
    if(Print_dir != "") {
      boost::filesystem::create_directory(Print_dir+"/extr");
      vector<double> X, Xerr, Y,Yerr, Extr_val;
      for(int imeas=0;imeas<Nmeas;imeas++) {
	X.push_back( Op[imeas].ave());
	Xerr.push_back( Op[imeas].err());
	Y.push_back( Meas[imeas].ave());
	Yerr.push_back( Meas[imeas].err());
	Extr_val.push_back(0);
      }
      
      X.push_back(Op_phys.ave());
      Xerr.push_back(Op_phys.err());
      Y.push_back(a0_fit.ave());
      Yerr.push_back(a0_fit.err());
      Extr_val.push_back(1);

      double der_val= (0.01*Op_phys*D_fit/a0_fit).ave();
      double der_err= (0.01*Op_phys*D_fit/a0_fit).err();
      
      Print_To_File({}, {X,Xerr, Y, Yerr, Extr_val}, Print_dir+"/extr/"+Tag+"_spline_extr.meas", "", "# X Xerr Y Yerr Extr_id, derivative: "+to_string_with_precision(der_val,6)+" +- "+to_string_with_precision(der_err,6));
      
      //print fitting function
      Vfloat ave_ops;
      for(auto & op: Op) ave_ops.push_back( op.ave());
      int Npoints = 300;
      double xmin= min(0.5*Op_phys.ave(), 0.5*(*min_element(ave_ops.begin(), ave_ops.end())));
      double xmax= max(1.5*Op_phys.ave(), 1.5*(*max_element(ave_ops.begin(), ave_ops.end())));
      Vfloat Xpoints;
      distr_t_list Fit_func;
      for(int ipoint=0;ipoint<Npoints;ipoint++) {
	double xp= xmin+ ipoint*(xmax-xmin)/((double)Npoints);
	if(Nmeas==2) { Fit_func.distr_list.push_back(a0_fit + D_fit*(xp-Op_phys)) ;}
	else if(Nmeas==3)  { Fit_func.distr_list.push_back(a0_fit + D_fit*(xp-Op_phys) + D2_fit*(xp-Op_phys.ave())*(xp-Op_phys));}
	else crash("In Meson_mass_extrapolation.cpp: SPLINATOR called with Nmeas: "+to_string(Nmeas));
	Xpoints.push_back(xp);
      }

    Print_To_File({}, {Xpoints, Fit_func.ave(), Fit_func.err()},Print_dir+"/extr/"+Tag+"_spline_extr.fit_func", "", "");
    


    }

    
  }

  return a0_fit;
  
}
