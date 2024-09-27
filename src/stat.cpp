#include "../include/stat.h"

using namespace std;


Pfloat JackAve(const Vfloat &JackDistr) {

  int Nclusters = JackDistr.size();
  Pfloat res= make_pair(0.0,0.0);
  for(auto & data: JackDistr) res.first += data/Nclusters;
  for(auto & data: JackDistr) res.second += pow((data - res.first),2);

  res.second = sqrt( ((Nclusters-1.0)/(double)Nclusters)*res.second);

 
  return res;
}

double Compute_jack_cov(const Vfloat& A,const Vfloat& B) {

  if(A.size() != B.size()) crash(" In Compute_jack_cov: Number of clusters for the two observables A,B don't match");

  int Nclusters = A.size();

  double barA = JackAve(A).first;
  double barB = JackAve(B).first;

  double result=0;

  for(int ijack=0;ijack<Nclusters;ijack++) result += (A[ijack]-barA)*(B[ijack]-barB);

  //Vfloat C = Multiply_vectors(A,B);
  //double AB = JackAve(C).first;
  return (Nclusters-1.0)*result/Nclusters;
  
}

void Compute_covariance_matrix(bool UseJack, Eigen::MatrixXd& Cov, int nargs,...) {

  
  Cov.resize(nargs,nargs);
  va_list args;
  va_start(args, nargs);
  VVfloat A(nargs);
  for(int i=0; i<nargs;i++) A[i] = va_arg(args,Vfloat);
  va_end(args);

  for(int i=0; i<nargs;i++)
    for(int j=i;j<nargs;j++) {
      double res;
      if(UseJack) res= Compute_jack_cov(A[i], A[j]);
      else res = Compute_boot_cov(A[i], A[j]);
      Cov(i,j) = res;
      Cov(j,i) = res;
    }

  return;
}


void Compute_correlation_matrix(bool UseJack, Eigen::MatrixXd& Corr, int nargs,...) {

  
  Corr.resize(nargs,nargs);
  va_list args;
  va_start(args, nargs);
  VVfloat A(nargs);
  for(int i=0; i<nargs;i++) A[i] = va_arg(args,Vfloat);
  va_end(args);

  for(int i=0; i<nargs;i++)
    for(int j=i;j<nargs;j++) {
      double num,den;
      if(UseJack) {
	num= Compute_jack_cov(A[i], A[j]);
	den = sqrt( Compute_jack_cov(A[i], A[i]) * Compute_jack_cov(A[j], A[j]));
      }
      else {
	num = Compute_boot_cov(A[i], A[j]);
	den = sqrt( Compute_boot_cov(A[i], A[i]) * Compute_boot_cov(A[j], A[j]));
      }
      Corr(i,j) = num/den;
      Corr(j,i) = num/den;
    }

  return;
}


void Compute_autocorrelation_time(const Vfloat &data, string Path, string Tag) {

  //compute empirical autocorrelation function

  int N = (signed)data.size();

  Vfloat rho_E(N, 0.0);
  Vfloat log_rho_E(N,0.0);
  Vfloat t_int(N,0.0);
  double t_accumulated = -1.0;

  for(int t=0; t< N;t++) {

    double bar_x_t=0;
    double bar_y_t=0;
    double num=0;
    double den=0;
    double den1=0;
    double den2=0;
    for(int i=0;i< N-t;i++) bar_x_t += (1.0/(double)(N-t))*data[i];
    for(int i=0;i< N-t;i++) bar_y_t += (1.0/(double)(N-t))*data[t+i];
    for(int i=0;i< N-t;i++) num += (data[i] - bar_x_t)*(data[t+i]-bar_y_t);
    for(int i=0;i< N-t;i++) {den1 += pow( data[i] - bar_x_t, 2); den2 += pow( data[t+i]-bar_y_t,2);}
    den = sqrt( den1*den2);
    rho_E[t] = num/den;
    log_rho_E[t] = log(fabs(num/den));
    t_accumulated += 2.0*num/den;
    t_int[t] = t_accumulated;


  }

  Print_To_File({}, {rho_E, log_rho_E, t_int}, Path+"/autocorr_"+Tag+".dat", "", "");



  return;
}

distr_t Jackknife::DoJack(function<double(const Vfloat&)> F, int Nobs,...) {


  distr_t JackDistribution(1);
  
  VVfloat Copy_data;

  va_list args;
  va_start(args, Nobs);
  for(int var_i=0;var_i<Nobs;var_i++) Copy_data.push_back(va_arg(args,Vfloat));
  va_end(args);

   //check that data size is same for all observables
  int N=0;
  for(unsigned int i=0; i<Copy_data.size();i++) {
    if((signed)Copy_data[i].size() != N && (signed)i != 0) crash("Jackknife::DoJack called for "+to_string(Nobs)+" observables having different data sizes");
    else N= Copy_data[i].size();
  }

  
  if(ThermalMeasToErase > 0) {
    for(auto &data: Copy_data) data.erase(data.begin(), data.begin()+ThermalMeasToErase);
  }

  if(Copy_data.size() == 0) crash("Jackknife routine called with an empty vector");


  if(Enable_fractional_jackknife) {

     //crash if called with Returnblocksizemax.
    if(ReturnBlockSizeMax) crash("DoJack called in mode fractional with Returnblocksizemax=true");
    double bs = ((double)N)/((double)Njacks); //fractional block_size
  
    //get total sum
    Vfloat total_sum(Nobs,0.0);
    for(int iobs=0;iobs<Nobs;iobs++)
      for(int iconf=0;iconf< N;iconf++) total_sum[iobs] += Copy_data[iobs][iconf];
    
    for(int ijack=0;ijack < Njacks; ijack++) {
      //get out of the block mean erasing the block ijack of size bs
      

      //initial bin time
      const double bin_start= ijack*bs;
      const double bin_end = bin_start + bs;

      //loop over time
      double binPos = bin_start;
      Vfloat data_to_cluster = total_sum;

      do {

	int iConf= floor(binPos + 1e-10);

	//Rectangle left point
	double lpoint = binPos;

	//Rectangle right point
	double rpoint = min(bin_end, iConf+1.0);

	//Rectangle horizontal size
	double rect_size = rpoint-lpoint;


	//add to jack
        for(int iobs=0; iobs <Nobs;iobs++) data_to_cluster[iobs] -= Copy_data[iobs][iConf]*rect_size;



	//update position
	binPos = rpoint;




      } while (bin_end - binPos > 1e-10);
      

      for(int iobs=0;iobs<Nobs;iobs++) data_to_cluster[iobs] /= (double)(N - bs);

      JackDistribution.distr.push_back(F(data_to_cluster));
    }
  }

  else {
    
  
  Vfloat sigma_obs_block(block_size_max,0);
  Vfloat mean_obs_block(block_size_max,0);
  Vfloat block_data;

  if(Verbose_jack) cout<<"#######JACKKNIFE ANALYSIS TEST MODE#######"<<endl;
  
  for(int block_size=1; block_size<=block_size_max;block_size++)
    {
      if(Verbose_jack) cout<<"block_size "<<block_size<<",block_size_max "<<block_size_max<<endl<<flush;
      if( (block_size==block_size_max && ReturnBlockSizeMax) || Verbose_jack || block_size == N/Njacks) {

	int num_bins;
	if( (block_size==block_size_max && ReturnBlockSizeMax) || Verbose_jack) num_bins= N/block_size;
	else num_bins = Njacks; //whether you want to fix the block size or the number of clusters to be used
	block_data.clear();
      
	for(int block_index=0;block_index<num_bins;block_index++) {


	  if(Verbose_jack) cout<<"block_index "<<block_index<<",num_bins "<<num_bins<<endl<<flush;
	  Vfloat data_to_cluster(Nobs,0);
	  for(int iconf=0; iconf<num_bins*block_size;iconf++ ) {
	    if( (iconf < block_index*block_size)  || ( iconf >= (block_index+1)*block_size) ) {
	      for(int iobs=0; iobs <Nobs;iobs++) data_to_cluster[iobs] += Copy_data[iobs][iconf]/(double)(num_bins*block_size-block_size);
	    }
	  }
	
	  
	  
	  block_data.push_back(F(data_to_cluster));
	  if(ReturnBlockSizeMax && block_size==block_size_max) JackDistribution.distr.push_back(F(data_to_cluster));
	  if(!ReturnBlockSizeMax && num_bins == Njacks) JackDistribution.distr.push_back(F(data_to_cluster));


	
	
	  
	}

	for(unsigned int scroll=0;scroll<block_data.size();scroll++)  mean_obs_block[block_size-1] += block_data[scroll];
	mean_obs_block[block_size-1] /= block_data.size();
     
      
	for(unsigned int scroll=0;scroll<block_data.size();scroll++) sigma_obs_block[block_size-1] += ((num_bins-1)/((double)num_bins))*pow((mean_obs_block[block_size-1] - block_data[scroll]),2);
	sigma_obs_block[block_size-1] = sqrt(sigma_obs_block[block_size-1]);

      
	if(Verbose_jack) cout<<block_size<<"\t"<<mean_obs_block[block_size-1]<<"\t"<<sigma_obs_block[block_size-1]<<endl<<flush;
      }
    }
  }

  return JackDistribution;
}





distr_t Jackknife::DoJack(int Nobs,...) {


  distr_t JackDistribution(1);
  
  VVfloat Copy_data;

  va_list args;
  va_start(args, Nobs);
  for(int var_i=0;var_i<Nobs;var_i++) Copy_data.push_back(va_arg(args,Vfloat));
  va_end(args);

  //check that data size is same for all observables
  int N=0;
  for(unsigned int i=0; i<Copy_data.size();i++) {
    if((signed)Copy_data[i].size() != N && (signed)i != 0) crash("Jackknife::DoJack called for "+to_string(Nobs)+" observables having different data sizes");
    else N= Copy_data[i].size();
  }

  
  if(ThermalMeasToErase > 0) {
    for(auto &data: Copy_data) data.erase(data.begin(), data.begin()+ThermalMeasToErase);
  }

  if(Copy_data.size() == 0) crash("Jackknife routine called with an empty vector");
  
  
  if(Enable_fractional_jackknife) {
    //crash if called with Returnblocksizemax.
    if(ReturnBlockSizeMax) crash("DoJack called in mode fractional with Returnblocksizemax=true");
    double bs = ((double)N)/((double)Njacks); //fractional block_size
 
    //get total sum
    double total_sum=0.0;
    for(int iobs=0;iobs<Nobs;iobs++)
      for(int iconf=0;iconf< N;iconf++) total_sum += Copy_data[iobs][iconf];
    
    for(int ijack=0;ijack < Njacks; ijack++) {
      //get out of the block mean erasing the block ijack of size bs


      //initial bin time
      const double bin_start= ijack*bs;
      const double bin_end = bin_start + bs;

      //loop over time
      double binPos = bin_start;
      double data_to_cluster = total_sum;

      do {

	int iConf= floor(binPos + 1e-10);

	//Rectangle left point
	double lpoint = binPos;

	//Rectangle right point
	double rpoint = min(bin_end, iConf+1.0);

	//Rectangle horizontal size
	double rect_size = rpoint-lpoint;


	//add to jack
        for(int iobs=0; iobs <Nobs;iobs++) data_to_cluster -= Copy_data[iobs][iConf]*rect_size;



	//update position
	binPos = rpoint;




      } while (bin_end - binPos > 1e-10);
      

      data_to_cluster /= (double)(N - bs);

      JackDistribution.distr.push_back(data_to_cluster);
    }


  }


  else {
  
  
  Vfloat sigma_obs_block(block_size_max,0);
  Vfloat mean_obs_block(block_size_max,0);
  Vfloat block_data;

  if(Verbose_jack) cout<<"#######JACKKNIFE ANALYSIS TEST MODE#######"<<endl;
  
  for(int block_size=1; block_size<=block_size_max;block_size++)
    {
      if(Verbose_jack) cout<<"block_size "<<block_size<<",block_size_max "<<block_size_max<<endl<<flush;
      if( (block_size==block_size_max && ReturnBlockSizeMax) || Verbose_jack || block_size == N/Njacks) {

	int num_bins;
	if( (block_size==block_size_max && ReturnBlockSizeMax) || Verbose_jack) num_bins= N/block_size;
	else num_bins = Njacks; //whether you want to fix the block size or the number of clusters to be used
	block_data.clear();
      
	for(int block_index=0;block_index<num_bins;block_index++) {


	  if(Verbose_jack) cout<<"block_index "<<block_index<<",num_bins "<<num_bins<<endl<<flush;
	  double data_to_cluster =0;
	  for(int iconf=0; iconf<num_bins*block_size;iconf++ ) {
	    if( (iconf < block_index*block_size)  || ( iconf >= (block_index+1)*block_size) ) {
	      for(int iobs=0; iobs <Nobs;iobs++) data_to_cluster += Copy_data[iobs][iconf];
	    }
	  }
	

	 


	  data_to_cluster /= (double)(num_bins*block_size-block_size); //out of the block mean
	  block_data.push_back(data_to_cluster);
	  if(ReturnBlockSizeMax && block_size==block_size_max) JackDistribution.distr.push_back(data_to_cluster);
	  if(!ReturnBlockSizeMax && num_bins == Njacks) JackDistribution.distr.push_back(data_to_cluster);


	
	
	  
	}

	for(unsigned int scroll=0;scroll<block_data.size();scroll++)  mean_obs_block[block_size-1] += block_data[scroll];
	mean_obs_block[block_size-1] /= block_data.size();
     
      
	for(unsigned int scroll=0;scroll<block_data.size();scroll++) sigma_obs_block[block_size-1] += ((num_bins-1)/((double)num_bins))*pow((mean_obs_block[block_size-1] - block_data[scroll]),2);
	sigma_obs_block[block_size-1] = sqrt(sigma_obs_block[block_size-1]);

      
	if(Verbose_jack) cout<<block_size<<"\t"<<mean_obs_block[block_size-1]<<"\t"<<sigma_obs_block[block_size-1]<<endl<<flush;
      }
    }

  }
  
  return JackDistribution;
}




distr_t Bootstrap::DoBoot(int Nobs, ...) {

  VVfloat Copy_data;
  va_list args;
  va_start(args,Nobs);
  for(int i=0;i<Nobs;i++) Copy_data.push_back(va_arg(args,Vfloat));
  va_end(args);
  distr_t boot_distribution(0, Nboots);
   
  for(int i=0; i<Nboots; i++) 
    for(unsigned int j=0; j<bootstrap_data[i].size();j++)
      for(int iobs=0;iobs<Nobs;iobs++) boot_distribution.distr[i] += Copy_data[iobs][bootstrap_data[i][j]]/Copy_data[iobs].size();
  
  
  return boot_distribution;  
}

distr_t Bootstrap::DoBoot(function<double(const Vfloat&)> F,int Nobs, ...) {

  VVfloat Copy_data;
  va_list args;
  va_start(args,Nobs);
  for(int i=0;i<Nobs;i++) Copy_data.push_back(va_arg(args,Vfloat));
  va_end(args);
  distr_t boot_distribution(0, Nboots);
   
  for(int i=0; i<Nboots; i++) {
    Vfloat Boot_obs(Nobs,0);
    for(unsigned int j=0; j<bootstrap_data[i].size();j++)
      for(int iobs=0;iobs<Nobs;iobs++) Boot_obs[iobs] += Copy_data[iobs][bootstrap_data[i][j]]/Copy_data[iobs].size();
    boot_distribution.distr[i] = F(Boot_obs);
  }
  
  return boot_distribution;  
}



Pfloat BootAve(const Vfloat &BootDistr) {

  Pfloat ave_err = make_pair(0.0,0.0);


  int N= BootDistr.size();
  for(int i=0; i<N; i++)  ave_err.first += BootDistr[i]/N;
  for(int i=0; i<N; i++) ave_err.second += pow( BootDistr[i]-ave_err.first,2)/(N-1.0);

  ave_err.second = sqrt(ave_err.second);
  return ave_err;
}



double Compute_boot_cov(const Vfloat& A,const Vfloat& B) {

    double barA = BootAve(A).first;
    double barB = BootAve(B).first;

  
    int N=A.size();
    double res=0.0;
    if(N != (signed)B.size()) crash("In Compute_boot_cov: A.size() != B.size()");
    for(int i=0;i<N;i++) res += (1.0/(N-1.0))*(A[i]-barA)*(B[i] - barB); 
    return res;
}




double distr_t::ave() const {
  
   if(UseJack) return JackAve(this->distr).first;

   return BootAve(distr).first;
 }
double distr_t::err() const {
  if(UseJack) {
    return JackAve(this->distr).second;
  }
  else return BootAve(this->distr).second;
}

Pfloat distr_t::ave_err() const {
   if(UseJack) return JackAve(this->distr);
   else return BootAve(distr);
 }
 
int distr_t::size() const {return distr.size();}

Vfloat distr_t_list::ave() const {
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(distr_list[i].ave());
    return res;
  }

double distr_t_list::ave(int i_distr) const {
    if(i_distr >= (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return distr_list[i_distr].ave();
  }
Vfloat distr_t_list::err()  const{
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].err());
    return res;
  }
double distr_t_list::err(int i_distr) const {
    if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return this->distr_list[i_distr].err();
  }

vector<Pfloat> distr_t_list::ave_err() const {
    vector<Pfloat> res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].ave_err());
    return res;
  }


  
Pfloat distr_t_list::ave_err(int i_distr) const{
     if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
     return distr_list[i_distr].ave_err();
  }

int distr_t_list::size() const { return distr_list.size();}

Vfloat distr_t_list::Get_distr_index(int k) const {
    Vfloat res;
    for(int i=0; i < this->size(); i++) {
      if(this->distr_list[i].size() <= k ) crash("In distr_t_list, call to I_list(int k) invalid. Positional argument k is too large");
      res.push_back(this->distr_list[i].distr[k]);
    }
    return res;
}

//operator overloading

distr_t operator*(const distr_t& A, const distr_t& B) {

  if(A.size() != B.size()) crash("In distr_t, call to A*B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B.distr[i];
   return res;
 }
 
 distr_t operator+(const distr_t& A,const  distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A+B is invalid. A and B have different sizes");
   distr_t res(A.UseJack, A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B.distr[i];
   return res;
 }
 
 distr_t operator-(const distr_t& A, const distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A-B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B.distr[i];
   return res;
 }

 distr_t operator/(const distr_t& A, const distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) {
     //if(fabs(B.distr[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element: "+to_string_with_precision(B.distr[i], 16));
     res.distr[i] = A.distr[i]/B.distr[i];
   }
   return res;
 }

 double operator%(const distr_t& A,const  distr_t& B) {

   if(A.UseJack != B.UseJack) crash("Error in distr_t%distr_t, the two distribution are not using the same method");
   if(A.UseJack) return Compute_jack_cov(A.distr, B.distr);
   else return Compute_boot_cov(A.distr, B.distr);
 }

  

 distr_t operator*(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B;
   return res;
 }

 distr_t operator*(double B,  const distr_t& A) {
   return A*B;
 }

 distr_t operator+(const distr_t& A, double B) {
  
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B;
   return res;
 }

 distr_t operator+(double B,const distr_t& A) {
   return A+B;
 }

 distr_t operator-(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B;
   return res;
 }

 distr_t operator-(double B,const distr_t& A) {
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = B-A.distr[i];
   return res;
 }

 distr_t operator/(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   //if(fabs(B)<eps(16)) crash("cannot call operator distr_t/double with double=0");
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]/B;
   return res;
 }

 distr_t operator/(double B,const distr_t& A) {
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) {
     //if(fabs(A.distr[i])<eps(16)) crash("Cannot call double/distr_t with distr_t containing zeros"); 
     res.distr[i] = B/A.distr[i];
   }
     return res;
 }

 distr_t operator*(const distr_t& A,const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A*B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B[i];
   return res;
 }
 
 distr_t operator+(const distr_t& A, const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A+B is invalid. A and B have different sizes");
   distr_t res(A.UseJack, A.size());
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B[i];
   return res;
 }
 
 distr_t operator-(const distr_t& A, const Vfloat& B) {
   
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A-B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B[i];
   return res;
 }
 
 distr_t operator/(const distr_t& A,const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) {
     //if(fabs(B[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element");
     res.distr[i] = A.distr[i]/B[i];
   }
   return res;
 }

 distr_t operator*(const Vfloat &B,const distr_t& A) {
   return A*B;
 }
 distr_t operator+(const Vfloat &B, const distr_t& A) {
   return A+B;
 }
 distr_t operator-(const Vfloat& B, const distr_t& A) {
   return -1.0*A + B;
 }
 distr_t operator/(const Vfloat& B, const distr_t& A) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) {
     //if(fabs(A.distr[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element");
     res.distr[i] = B[i]/A.distr[i];
   }
   return res;
 }


  
distr_t_list operator*(const distr_t_list& A, const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A*B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B.distr_list[i];
    return res;
}

distr_t_list operator+(const distr_t_list& A,const distr_t_list& B) {

  if(A.size() != B.size()) crash("In distr_t_list, call to A+B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B.distr_list[i];
  return res;
}

distr_t_list operator-(const distr_t_list& A, const distr_t_list& B) {

  if(A.size() != B.size()) crash("In distr_t_list, call to A-B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A,const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A/B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B.distr_list[i];
    
  return res;
}

Vfloat operator%(const distr_t_list& A,const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A%B is invalid. A and B have different sizes");
  Vfloat cov(A.size());
  for(int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.distr_list[i].size()) crash("In distr_t_list, call to A%B invalid, the two operands have different distr_t sizes");
    cov[i] = A.distr_list[i]%B.distr_list[i];
  }
  return cov;
}

////////////////////////////////

distr_t_list operator*(const distr_t_list& A, const distr_t& B) {

  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B;
  return res;
}
distr_t_list operator*(const distr_t& B, const distr_t_list& A) {return A*B;}

distr_t_list operator+(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
  return res;
}
distr_t_list operator+(const distr_t& B, const distr_t_list& A) { return A+B;}

distr_t_list operator-(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
}

distr_t_list operator-(const distr_t& B,const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}

distr_t_list operator/(const distr_t& B,const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}
  

Vfloat operator%(const distr_t_list& A,const distr_t& B) {
  Vfloat cov(A.size());
  for( int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.size()) crash("In distr_t_list, call to A%B<distr_t_list, distr_t>  invalid, the two operands have different distr_t sizes");
    cov[i] = A.distr_list[i]%B;
  }
  return cov;
}

Vfloat operator%(const distr_t& B, const distr_t_list& A) {
  return A%B;
}

///////////////////////////////


distr_t_list operator*(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
    for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]*B;
    return res;
}

distr_t_list operator*(double B, const distr_t_list& A) {return A*B;}

distr_t_list operator+(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
     return res;
}

distr_t_list operator+(double B, const distr_t_list& A) { return A+B;}

distr_t_list operator-(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
  }

distr_t_list operator-(double B, const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}
  
distr_t_list operator/(double B, const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}

distr_t_list operator+(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]+B[i];

  return res;
}

distr_t_list operator-(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]-B[i];

  return res;
}

distr_t_list operator+(const Vfloat &B, const distr_t_list &A) { return A + B; }
distr_t_list operator-(const Vfloat& B, const distr_t_list& A) { return -1.0*(A-B); }


distr_t_list operator*(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]*B[i];

  return res;
}

distr_t_list operator*(const Vfloat& B, const distr_t_list& A) { return A*B;}


distr_t Get_id_jack_distr(int N) {


  distr_t id_jack_distr;

  for(int i=0;i<N;i++) id_jack_distr.distr.push_back(1);

  return id_jack_distr;
}

distr_t Get_id_distr(int N, bool UseJack) {


  distr_t id_jack_distr(UseJack);

  for(int i=0;i<N;i++) id_jack_distr.distr.push_back(1);

  return id_jack_distr;
}

distr_t_list Get_id_jack_distr_list(int size, int N) {
  distr_t_list return_distr(1);
  for(int i=0; i<size;i++) return_distr.distr_list.push_back( Get_id_jack_distr(N));
  return return_distr;
}

distr_t_list Get_id_distr_list(int size, int N, bool UseJack) {
  distr_t_list return_distr(UseJack);
  for(int i=0; i<size;i++) return_distr.distr_list.push_back( Get_id_distr(N,UseJack));
  return return_distr;
}


distr_t AIC( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas,  bool NEIL, const vector<double> &mult ) {

  Vfloat MULT;
  if(mult.size() == 1 && mult[0] == 0) { MULT.resize(VAL.size(),1);}
  else MULT=mult;
  
  if( (VAL.size() != ch2.size()) || (VAL.size() != Ndof.size()) || ( Ndof.size() != Nmeas.size())) crash("AIC analysis called with vectors of different sizes");
  
  int N=VAL.size();
  assert(N > 0);
  double F= ((NEIL)?2:1);

  
  distr_t RES= 0.0*VAL[0];
  
  vector<double> W(N,0.0);
  double wtot=0; double syst=0;
  for(int i=0; i<N;i++) { W[i] = MULT[i]*exp(-0.5*( ch2[i] + 2*( Nmeas[i] - Ndof[i]) - F*Nmeas[i])); wtot += W[i]; }
  for(int i=0; i<N;i++) { W[i] /= wtot; RES = RES + W[i]*VAL[i];}
  for(int i=0; i<N;i++) {syst += W[i]*pow( VAL[i].ave() - RES.ave(),2); }
  syst= sqrt(syst);

  RES = RES.ave() + (sqrt(pow(syst,2)+pow(RES.err(),2))/RES.err())*(RES - RES.ave());
 
 
  return RES;
}


distr_t AIC_lin( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas,  bool NEIL, const vector<double> &mult ) {

  Vfloat MULT;
  if(mult.size() == 1 && mult[0] == 0) { MULT.resize(VAL.size(),1);}
  else MULT=mult;

    
  if( (VAL.size() != ch2.size()) || (VAL.size() != Ndof.size()) || ( Ndof.size() != Nmeas.size())) crash("AIC analysis called with vectors of different sizes");

  int N=VAL.size();
  assert(N > 0);
  double F= ((NEIL)?2:1);
   
  distr_t RES= 0.0*VAL[0];

  vector<double> W(N,0.0);
  double wtot=0; double syst=0;
  for(int i=0; i<N;i++) { W[i] = MULT[i]*exp(-0.5*( ch2[i] + 2*( Nmeas[i] - Ndof[i]) - F*Nmeas[i])); wtot += W[i]; }
  for(int i=0; i<N;i++) { W[i] /= wtot; RES = RES + W[i]*VAL[i];}
  for(int i=0; i<N;i++) {syst += W[i]*pow( VAL[i].ave() - RES.ave(),2); }
  syst= sqrt(syst);

  RES = RES + (syst/RES.err())*(RES - RES.ave());
 
 
  return RES;
  
   

}

distr_t BMA_uniform( const vector<distr_t> &VAL, const vector<double> &ch2, const vector<int> &Ndof, const vector<int> &Nmeas, double ch2_th) {


  
  if( (VAL.size() != ch2.size()) || (VAL.size() != Ndof.size()) || ( Ndof.size() != Nmeas.size())) crash("AIC analysis called with vectors of different sizes");

  int N=VAL.size();
  
  distr_t RES= 0.0*VAL[0];
  
  vector<double> W(N,0.0);
  double wtot=0; double syst=0;
  for(int i=0; i<N;i++) { W[i] = ((ch2[i]/Ndof[i] < ch2_th)?1.0:0) ; wtot += W[i]; }
  for(int i=0; i<N;i++) { W[i] /= wtot; RES = RES + W[i]*VAL[i];}
  for(int i=0; i<N;i++) {syst += W[i]*pow( VAL[i].ave() - RES.ave(),2); }
  syst= sqrt(syst);

  RES = RES.ave() + (sqrt(pow(syst,2)+pow(RES.err(),2))/RES.err())*(RES - RES.ave());
 
 
  return RES;
}

distr_t BMA_Eq_29( const vector<distr_t> &VAL) {


  
 
  int N=VAL.size();
  
  distr_t RES= 0.0*VAL[0];
  
  vector<double> W(N,0.0);
  double wtot=0; double syst=0;
  for(int i=0; i<N;i++) { W[i] = 1.0/N ; wtot += W[i]; }
  for(int i=0; i<N;i++) { W[i] /= wtot; RES = RES + W[i]*VAL[i];}
  for(int i=0; i<N;i++) {syst += W[i]*pow( VAL[i].ave() - RES.ave(),2); }
  syst= sqrt(syst);

  RES = RES.ave() + (sqrt(pow(syst,2)+pow(RES.err(),2))/RES.err())*(RES - RES.ave());
 
 
  return RES;


}




distr_t_list EXP_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list EXPT_DL(const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr_list[t].distr[i]*t);
  return B;
}


distr_t_list COSH_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = cosh(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list SINH_DL (const distr_t_list & A) { distr_t_list B=A;
      for(int t=0;t<A.size();t++)
	for(int i=0;i<A.distr_list[t].size();i++) B.distr_list[t].distr[i] = sinh(A.distr_list[t].distr[i]);
      return B;
 }

distr_t_list SQRT_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = sqrt(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list LOG_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = log(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list POW_DL (const distr_t_list &A, int n) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = pow(A.distr_list[t].distr[i],n);
  return B;
}


distr_t_list EXPT_D(const distr_t &A, int T) {
  distr_t_list B(A.UseJack, T, A.size()); 
  for (int t = 0; t < T; t++)
    for (int i = 0; i < A.size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr[i]*t);
  return B;
}


distr_t EXP_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = exp(A.distr[i]);
  return B;
}

distr_t COSH_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = cosh(A.distr[i]);
  return B;
}

distr_t SINH_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = sinh(A.distr[i]);
  return B;
}

distr_t LOG_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = log(A.distr[i]);
  return B;
}

distr_t SQRT_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = sqrt(A.distr[i]);
  return B;
}

distr_t POW_D (const distr_t  &A,int n) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = pow(A.distr[i],n);
  return B;
}


