#ifndef __input__
#define __input__

#include "numerics.h"


using namespace std;


class data_t {

 public:
  data_t() {Directory=""; Obs_name=""; File_name=""; SC="#";}
  data_t(string A, string B, string C) : Directory(A), File_name(B), Obs_name(C) {SC ="#"; }
  data_t(string A, string B, string C, string D) : Directory(A), File_name(B), Obs_name(C), SC(D) {}

  void Read();
  void Read(string A, string B, string C) { Directory=A; File_name=B; Obs_name=C; Read();}
  void Read(string A, string B, string C, string D) { Directory=A; File_name=B; Obs_name=C; SC=D; Read(); }
  VVVfloat col(int icol) {
    map<int, VVVfloat>::iterator it;
    it= data.find(icol);
    if(it == data.end()) crash("In class data_t, column: "+to_string(icol)+" is indefined");
    return it->second;
  }
  

  vector<string> Tag;
  map<int, VVVfloat> data; //corr<col, [ens][t][conf]
  int size;
  Vint Nconfs;
  Vint nrows;
  int Ncols;
  string Directory;
  string File_name;
  string Obs_name;
  string SC;
};

class file_t {

 public:

 file_t() : nrows(0), path(""), SC("#"), raw_data(0) {} 
 file_t(string file_path) : nrows(0), path(file_path), SC("#"), raw_data(0) {}


  void Read();
  void Read(string file_path) { path=file_path; Read();};

  Vfloat col(int icol) { return data[icol];}


  double nrows;
  string path;
  string SC;
  vector<string> raw_data;
  VVfloat data;
 



};








#endif
