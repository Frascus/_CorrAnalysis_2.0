#include "../include/pt3_momenta.h"

using namespace std;


void Add_to_mom_list(pt3_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt3_momenta> A;


  for(auto & c: header.comb) A.emplace_back(Vfloat{c.th0[0], c.th0[1], c.th0[2]},Vfloat{c.ths[0], c.ths[1], c.ths[2]},Vfloat{c.tht[0], c.tht[1], c.tht[2] }, c.mu1, c.mu2, c.off,L, header.tmax); 
  M.mom.push_back(A);
  return;
}


void Add_to_mom_list(pt2_momenta_list &M, struct header_virph &header, double& L) {

  vector<pt2_momenta> A;

  for(auto & c: header.comb) A.emplace_back( Vfloat{c.th0[0], c.th0[1], c.th0[2]}, Vfloat{c.ths[0], c.ths[1], c.ths[2]}, c.mu1, c.mu2, L, c.i0, c.is);

  M.mom.push_back(A);

  return;
}
