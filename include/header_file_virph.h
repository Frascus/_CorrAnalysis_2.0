#ifndef __header_file_virph__
#define __header_file_virph__

#include "numerics.h"
#include "binary_io.h"

using namespace std;

struct flavour_t
{
  int qhat;
  double kappa;
  double mu;
  double su3csw;
  double u1csw;
  double cF;
  double cF_prime;
  double th1;
  double th2;
  double th3;
};

struct inv_t
{
  int isolv;
  double mu,th[3];
};

struct einv_t
{
  int isolv;
  double mu,th0[3],tht[3],off;
};

struct combination_t
{
  int i0,it,is;
  double mu1,mu2,off;
  double th0[3],tht[3],ths[3];
};

struct  header_virph
{
  int twist;
  int nf;
  int nsrc,nsrcd;
  int l0,l1,l2,l3;
  int nk,nmoms;
  double beta,ksea,musea,csw;
  double *k,*mu,**mom;
  int allocated=0; 
  int tmax;
  int x0;
  int stype;
  int ninv;
  int neinv;
  int nsolv;
  int nhits;
  int phptype;
  int z0;
  int ncomb;
  int ngsm;
  double epsgsm;
  int nqsml,nqsm0,nqsm;
  double epsqsm;
  flavour_t gflv;
  vector<combination_t> comb;
  vector<inv_t> inv;
  vector<einv_t> einv;
  
  
  int header_size;
  int file_size;
  int file_nconf;
} ;


void read_header_bin(FILE *stream, struct header_virph &header);

int Get_symmetric_comb(struct header_virph &header, int icomb) ;

int Get_comb_k0(struct  header_virph &header, int icomb);

int Get_comb_k0_same_off(struct header_virph &header, int icomb);

int Get_2pt_k0p0(struct header_virph &header,double mu1, double mu2);

int Get_2pt_p(struct header_virph &header, int i0, int is);

VVfloat Get_obs_2pt(FILE *stream, struct header_virph &header, int ire,int icomb, int icorr, int smearing_level);


VVfloat Get_obs_3pt(FILE *stream, struct header_virph &header, int ire, int icomb, int alpha, int mu, string A, int smearing_level);

int Get_number_of_configs_3pt(FILE* stream, struct header_virph& header);

int Get_number_of_configs_2pt(FILE* stream, struct header_virph& header);


#endif
