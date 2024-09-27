#ifndef __input__
#define __input__

#include "numerics.h"


using namespace std;


class data_t {

 public:
  data_t() {Directory=""; Obs_name=""; File_name=""; SC="#"; sort_Custom=false; CustomSorting= [](string A, string B) -> bool { return A<B;};}
  data_t(string A, string B, string C) : Directory(A), File_name(B), Obs_name(C) {SC ="#";  sort_Custom=false; CustomSorting= [](string A, string B) -> bool { return A<B;};}
  data_t(string A, string B, string C, string D) : Directory(A), File_name(B), Obs_name(C), SC(D) {sort_Custom=false; CustomSorting= [](string A, string B) -> bool { return A<B;};}
  data_t(string A, string B, string C,const  function<bool(string A, string B)> &F)  : Directory(A), File_name(B), Obs_name(C) {SC="#"; sort_Custom=true; CustomSorting=F;} 
  data_t(string A, string B, string C, string D, const function<bool(string A, string B)> &F)  : Directory(A), File_name(B), Obs_name(C), SC(D) { sort_Custom=true; CustomSorting=F;} 

  void Read();
  void Read(string A, string B, string C) { Directory=A; File_name=B; Obs_name=C; Read();}
  void Read(string A, string B, string C, string D) { Directory=A; File_name=B; Obs_name=C; SC=D; Read(); }
  void Read(string A, string B, string C,const function<bool(string A, string B)> &F) {Directory=A; File_name=B; Obs_name=C; sort_Custom=true; CustomSorting=F; Read();}
  void Read(string A, string B, string C, string D,const  function<bool(string A, string B)> &F) {Directory=A; File_name=B; Obs_name=C;SC=D; sort_Custom=true; CustomSorting=F; Read();}
  VVVfloat col(int icol) {
    map<int, VVVfloat>::iterator it;
    it= data.find(icol);
    if(it == data.end()) crash("In class data_t, column: "+to_string(icol)+" is indefined");
    return it->second;
  }
  int Get_iens_from_tag(string T) {
    for(unsigned int iens=0;iens< Tag.size();iens++) {
      if(T==Tag[iens]) return iens;
    }
    crash("Cannot find Ensemble Tag: "+T+" in Tag list");
    return 0;
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
  bool sort_Custom;
  function<bool(string A,string B)> CustomSorting;
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


Vfloat Read_From_File(string Path, int icol, int ncols, int rows_to_skip=0);




#endif
