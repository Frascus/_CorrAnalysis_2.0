#include "../include/header_file_virph.h"


using namespace std;

const double e= 1.0; 


void read_header_bin(FILE *stream, struct header_virph &header)
{

   
    for(auto& i : {&header.tmax,&header.x0,&header.stype,&header.phptype,&header.z0,&header.ninv,&header.neinv,&header.nsolv,&header.nhits,&header.ncomb,&header.ngsm,&header.nqsml, &header.nqsm0,&header.nqsm})
        bin_read(*i,stream);


 
    header.inv.resize(header.ninv);
    header.comb.resize(header.ncomb);
    header.einv.resize(header.neinv);

    fread(&header.epsgsm,sizeof(double),1,stream);
    fread(&header.epsqsm,sizeof(double),1,stream);
    
    fread(&header.gflv.qhat,sizeof(int),1,stream);


    
    for(auto& d : {&header.gflv.kappa,&header.gflv.mu,&header.gflv.su3csw,&header.gflv.u1csw,&header.gflv.cF,&header.gflv.cF_prime,&header.gflv.th1,&header.gflv.th2,&header.gflv.th3})
        bin_read(*d,stream);


    
    for(int icomb=0;icomb<header.ncomb;++icomb)
        {
        auto& c=header.comb[icomb];
        for(auto& i : {&c.i0,&c.it,&c.is})
            bin_read(*i,stream);
        for(auto& d : {&c.mu1,&c.mu2,&c.off})
            bin_read(*d,stream);
        for(auto& d : {&c.th0[0],&c.th0[1],&c.th0[2]})
            bin_read(*d,stream);
        for(auto& d : {&c.tht[0],&c.tht[1],&c.tht[2]})
            bin_read(*d,stream);
        for(auto& d : {&c.ths[0],&c.ths[1],&c.ths[2]})
            bin_read(*d,stream);
    }
    
    for(int inv=0;inv<header.ninv;++inv)
    {
      auto& v=header.inv[inv];
      for(auto& d : {&v.mu,&v.th[0],&v.th[1],&v.th[2]})
            bin_read(*d,stream);
    }
    header.header_size=ftell(stream);      
        
  
}







int Get_symmetric_comb(struct header_virph &header, int icomb) {

    int ci=-1;
    int i0=header.comb[icomb].i0;
    int is=header.comb[icomb].is;
    int found=0;
    
     
    for (int i =0; i<header.ncomb;i++ ){
        if (header.comb[i].i0==is   &&  header.comb[i].is==i0  )
        if (fabs(header.comb[i].off-header.comb[icomb].off )<1e-10)
        if (fabs(header.comb[i].th0[0]-header.comb[icomb].th0[0] )<1e-10)    
        if (fabs(header.comb[i].tht[0]-header.comb[icomb].tht[0] )<1e-10)        
        if (fabs(header.comb[i].ths[0]-header.comb[icomb].ths[0] )<1e-10)
        if (fabs(header.comb[i].th0[1]-header.comb[icomb].th0[1] )<1e-10)    
        if (fabs(header.comb[i].tht[1]-header.comb[icomb].tht[1] )<1e-10)        
        if (fabs(header.comb[i].ths[1]-header.comb[icomb].ths[1] )<1e-10)
        if (fabs(header.comb[i].th0[2]-header.comb[icomb].th0[2] )<1e-10)    
        if (fabs(header.comb[i].tht[2]-header.comb[icomb].tht[2] )<1e-10)        
        if (fabs(header.comb[i].ths[2]-header.comb[icomb].ths[2] )<1e-10)
	if (fabs(header.comb[i].off - header.comb[icomb].off ) < 1e-10)  {
                    ci=i;
                    found++;
        }
    }

    
			  
    if (found!=1){
    auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_opposite_mu\n");
	printf("ci: %d\n", ci);
        printf("Either there is no combination with opposite mu either there are many\n");
        exit(3);
    }

 
    return ci;

}




int Get_comb_k0(struct  header_virph &header, int icomb){
    int ci=-1;
    int found=0;
    int foundk;
          
    for (int i =0; i<header.ncomb;i++ ){
        foundk=0;
        for (int k=0; k<3;k++){
            if (fabs(header.comb[i].th0[k]-header.comb[i].tht[k] )<1e-7){
                if (fabs(header.comb[i].th0[k]-header.comb[icomb].th0[k] )<1e-7)
                if (fabs(header.comb[i].ths[k]-header.comb[icomb].ths[k] )<1e-7)
                if (fabs(header.comb[i].mu1-header.comb[icomb].mu1 )<1e-7)    
                if (fabs(header.comb[i].mu2-header.comb[icomb].mu2 )<1e-7)        
                if (fabs(header.comb[i].off)<1e-8){
                    foundk++;
                }
            }
                    
        }
        if (foundk==3){
            ci=i;
            found++;
        }
          
    }
    
    if (found!=1){
        auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_k0\n");
	printf("ci: %d\n", ci);
        printf("Either there is no combination with opposite mu either there are many\n");
        exit(3);
    }
   
        
    return ci;
}

int Get_comb_k0_same_off(struct  header_virph &header, int icomb){
    int ci=-1;
    int found=0;
    int foundk;
          
    for (int i =0; i<header.ncomb;i++ ){
        foundk=0;
        for (int k=0; k<3;k++){
            if (fabs(header.comb[i].th0[k]-header.comb[i].tht[k] )<1e-7){
                if (fabs(header.comb[i].th0[k]-header.comb[icomb].th0[k] )<1e-7)
                if (fabs(header.comb[i].ths[k]-header.comb[icomb].ths[k] )<1e-7)
                if (fabs(header.comb[i].mu1-header.comb[icomb].mu1 )<1e-7)    
                if (fabs(header.comb[i].mu2-header.comb[icomb].mu2 )<1e-7)        
		  if (fabs(header.comb[i].off- header.comb[icomb].off) < 1e-7){
                    foundk++;
                }
            }
                    
        }
        if (foundk==3){
            ci=i;
            found++;
        }
          
    }
    
    if (found!=1){
        auto c=header.comb[icomb];
        printf("icomb=%d\n",icomb);
        printf("i0 it is=%d  %d  %d\n",c.i0,c.it,c.is);
        printf("mu0  mut  off=%f  %f  %f\n",c.mu1,c.mu2,c.off);
        printf("th0=%f  %f  %f\n",c.th0[0],c.th0[1],c.th0[2]);
        printf("tht=%f  %f  %f\n",c.tht[0],c.tht[1],c.tht[2]);
        printf("ths=%f  %f  %f\n",c.ths[0],c.ths[1],c.ths[2]);
        printf("find_icomb_with_k0\n");
	printf("ci: %d\n", ci);
        printf("Get_comb_k0_same_off: Error in find the combination with kz=0 and same offshellness");
    }

    if( found != 1) return -1;
        
    return ci;
} 

int Get_2pt_k0p0(struct header_virph &header,double mu1, double mu2) {

  for(int icomb=0; icomb<header.ncomb; icomb++) {
    auto c = header.comb[icomb];

    if( fabs(c.mu1 -mu1) < eps(15) && fabs(c.mu2 -mu2)< eps(15)) 
      if( (pow(c.th0[0],2) + pow(c.th0[1],2) + pow(c.th0[2],2) + pow(c.ths[0],2) + pow(c.ths[1],2) + pow(c.ths[2],2) < 1.0e-7) && (c.off < 1.0e-8)) return icomb;

  }

  crash("In Get_2pt_k0p0 cannot find the 2pt function with zero momentum");
  return 0;
}


int Get_2pt_p(struct header_virph &header, int i0, int is) {


  for(int icomb=0; icomb< header.ncomb; icomb++) {
    auto c = header.comb[icomb];
    if (c.is==is && c.i0==i0) return icomb;
  }

  crash("In Get_2pt_p cannot find the 2pt function with opposite (i0, is) ");
  return 0;
}



VVfloat Get_obs_2pt(FILE *stream, struct header_virph &header, int ire,int icomb, int icorr, int smearing_level) {


  int tmp= header.header_size;

  int ncorr=5;

  
  int Nt = header.tmax;

 

  if(smearing_level >= header.nqsml) crash("In Get_obs_2pt: smearing level: "+to_string(smearing_level)+" not yet present. Aborting.");

  if(icomb >= header.ncomb) crash("In Get_obs_2pt: icomb > ncomb");


  auto c = header.comb[icomb];


  fseek(stream, 0, SEEK_END);
  int file_size_bytes = ftell(stream);
  int block_size_bytes = sizeof(int) + sizeof(double)*(2*header.tmax*header.ninv*header.ninv*header.nqsml*ncorr);
  int nconfs = (file_size_bytes-tmp)/block_size_bytes;
  if( (file_size_bytes-tmp)%block_size_bytes != 0) {
    cout<<"file size bytes: "<<file_size_bytes<<endl;
    cout<<"block size bytes: "<<block_size_bytes<<endl;
    cout<<"header size: "<<tmp<<endl;
    cout<<"nconfs: "<<nconfs<<endl;
    crash(" In Get_obs_2pt block_data_size does not divide file_size- header_size");
  }
  
  VVfloat obs(Nt);
  for(auto &o: obs) o.resize(nconfs);


  for(int iconf=0; iconf<nconfs;iconf++) {
    for(int t=0; t< header.tmax; t++) {
      double read=0;
      int pos = sizeof(int) +iconf*block_size_bytes;
      int loc_pos = sizeof(double)*(ire + 2*(t + Nt*(smearing_level + header.nqsml*(c.i0 + header.ninv*(c.is + header.ninv*icorr)))));

      fseek(stream, tmp+pos+loc_pos, SEEK_SET);
      bin_read(read, stream);
      obs[t][iconf] = read;
    }
  }


  return obs;   
      
  
}
  


VVfloat Get_obs_3pt(FILE *stream, struct header_virph &header, int ire, int icomb, int alpha, int mu, string A, int smearing_level) {

  int icorr;
  if(A=="V") icorr=1;
  else if(A=="A") icorr=0;
  else crash("In Get_obs_3pt, asking for a correlator that is neither Axial or Vector. Aborting.");

  int ncorr=2;

  int mus=4;
  int alphas=4;
  
  int ndim=4;

  int Nt= header.tmax;

  if(smearing_level >= header.nqsml) crash("In Get_obs_3pt: smearing level: "+to_string(smearing_level)+" not yet present. Aborting. ");

  if(icomb >= header.ncomb) crash("In Get_obs_3pt: icomb > ncomb");

  fseek(stream, 0, SEEK_END);
  int file_size_bytes = ftell(stream);
  int tmp= header.header_size;
  int block_size_bytes = sizeof(int) + sizeof(double)*(2*header.tmax*alphas*mus*header.ncomb*header.nqsml*ncorr);
  int nconfs = (file_size_bytes-tmp)/block_size_bytes;
 
  if( (file_size_bytes-tmp)%block_size_bytes != 0) crash(" In Get_obs_3pt block_data_size does not divide file_size- header_size");

  VVfloat obs(Nt);

  for(auto &o: obs) o.resize(nconfs);


  for(int t=0; t< Nt; t++) {  
    for(int iconf=0; iconf<nconfs;iconf++) {
      
      double read;
      int pos = sizeof(int) + iconf*block_size_bytes;
      int loc_pos = sizeof(double)*(ire + 2*(t + Nt*(alpha + ndim*(mu + ndim*( smearing_level + header.nqsml*(icomb + header.ncomb*icorr))))));

      fseek(stream, tmp+pos+loc_pos, SEEK_SET);
      bin_read(read, stream);
      obs[t][iconf] = e*read;
    }
  }
  

 return obs;
 
}


int Get_number_of_configs_3pt(FILE* stream, struct header_virph& header) {

 int ncorr=2;

 int mus=4;
 int alphas=4;
  
 fseek(stream, 0, SEEK_END);
 int file_size_bytes = ftell(stream);
 int tmp= header.header_size;
 int block_size_bytes = sizeof(int) + sizeof(double)*(2*header.tmax*alphas*mus*header.ncomb*header.nqsml*ncorr);
 int nconfs = (file_size_bytes-tmp)/block_size_bytes;
 
  if( (file_size_bytes-tmp)%block_size_bytes != 0) crash(" In Get_number_of_configs_3pt block_data_size does not divide file_size- header_size");

  return nconfs;

}


int Get_number_of_configs_2pt(FILE* stream, struct header_virph& header) {
  
  int tmp= header.header_size;

  int ncorr=5;

  fseek(stream, 0, SEEK_END);
  int file_size_bytes = ftell(stream);
  int block_size_bytes = sizeof(int) + sizeof(double)*(2*header.tmax*header.ninv*header.ninv*header.nqsml*ncorr);
  int nconfs = (file_size_bytes-tmp)/block_size_bytes;
  if( (file_size_bytes-tmp)%block_size_bytes != 0) {
    cout<<endl;
    cout<<"header.tmax: "<<header.tmax<<endl;
    cout<<"header.ninv: "<<header.ninv<<endl;
    cout<<"header.neinv: "<<header.neinv<<endl;
    cout<<"header.ncomb: "<<header.ncomb<<endl;
    cout<<"header.nqsml: "<<header.nqsml<<endl;
    cout<<"file size bytes: "<<file_size_bytes<<endl;
    cout<<"block size bytes: "<<block_size_bytes<<endl;
    cout<<"header size: "<<tmp<<endl;
    cout<<"nconfs: "<<nconfs<<endl;
    crash(" In Get_number_of_configs_2pt block_data_size does not divide file_size- header_size");
  }
    
  return nconfs;
}
