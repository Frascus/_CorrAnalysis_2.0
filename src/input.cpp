#include "../include/input.h"


using namespace std;



void data_t::Read() {

  Tag.clear();
  Nconfs.clear();
  nrows.clear();
  data.clear();
  

  vector<vector<vector<string>>> raw_data;


  //reading ensemble list
  boost::filesystem::directory_iterator end_itr;

  for (boost::filesystem::directory_iterator itr(Directory); itr != end_itr; ++itr) {
  

    if(!boost::filesystem::is_regular_file(itr->path())) {

      Tag.push_back(itr->path().string().substr(Directory.length()+1));
    }
  }

  size= Tag.size();


  //loop over ensembles
  int it=0, ncols_old=0; //sanity checks
  for(auto & ens: Tag) {
    int iconf=0;
    vector<string> Confs;
    vector<vector<string>> raw_data_ens;
    for(boost::filesystem::directory_iterator itr(Directory+"/"+ens); itr != end_itr; ++itr) {
      if(!boost::filesystem::is_regular_file(itr->path())) {Confs.push_back(itr->path().string()); iconf++;}
    }
    Nconfs.push_back(iconf);
   
    if(!sort_Custom) sort(Confs.begin(), Confs.end());
    else sort(Confs.begin(), Confs.end(), CustomSorting);

    //test
    //cout<<"###########"<<endl;
    //for(int ht=0;ht<(signed)Confs.size();ht++) cout<<Confs[ht]<<endl;
    //
    //cout<<"###########"<<endl;
    //
    
    int N_rows_old=0;
    int N_rows=0;
    int it2=0;
    for(auto &conf_path: Confs) {
      N_rows=0;
      vector<string> raw_data_conf;
      if(!boost::filesystem::is_regular_file(conf_path+"/"+File_name)) crash("Cannot open file_name: "+conf_path+"/"+File_name);
	 ifstream Read(conf_path+"/"+File_name);
	 string ReadLine="";
	 if(Obs_name != "") {while(ReadLine != Obs_name) {Read>>ReadLine; if(Read.eof()) crash("Could not find observable: "+Obs_name+" in file: "+File_name);}}
	 //read file
	 do { 
	   do {getline(Read, ReadLine);} while(ReadLine=="" && !Read.eof());
	   if(!Read.eof() && ReadLine.find(SC) == string::npos) {raw_data_conf.push_back(ReadLine); N_rows++;}
	   int ncols=0;
	   stringstream s(ReadLine);
	   double tmp=0;
	   while(s>>tmp) ncols++;
	   if(ncols != ncols_old && it!=0 && !Read.eof() && (ReadLine.find(SC) == string::npos)) {
	     cout<<"Number of columns in file: "+File_name+", obs: "+Obs_name+" is not constant over configs"<<endl;
	     cout<<"ncols_old: "<<ncols_old<<" ncols: "<<ncols<<endl;
	     cout<<"Ens tag: "<<ens<<", conf_path: "<<conf_path<<endl;
	     cout<<"iteration: "<<it<<endl;
	     cout<<"Conf num: "<<it2<<endl;
	     crash("Exiting...");
	   }
	   if(!Read.eof() && (ReadLine.find(SC) == string::npos)) {
	     Ncols=ncols;
	     ncols_old=ncols;
	     it++;
	   }
	 } while( !Read.eof() && (ReadLine.find(SC) == string::npos || Obs_name == ""));

           
	 if(N_rows != N_rows_old && it2!=0) {
	   cout<<"While reading obs: "+Obs_name+" for ensemble: "+ens+" number of rows between configurations do not match "<<endl;
	   cout<<"N_rows_old= "<<N_rows_old<<" N_rows: "<<N_rows<<endl;
	   crash("Exiting...");
	 }
	 N_rows_old = N_rows;
	 it2++;
	 raw_data_ens.push_back(raw_data_conf);
	 }
      nrows.push_back(N_rows);
      raw_data.push_back(raw_data_ens);
  }
  

     
  
    //split the string into several double and assign them  to R.data
    vector<VVVfloat> convert_raw_data(Ncols);
    for(auto &c1: convert_raw_data) {
      c1.resize(this->size);
      for(int i=0; i<this->size;i++)  {
	c1[i].resize(Nconfs[i]);
      }
    }
  

    for(int i=0; i<this->size;i++)
      for(int j=0; j<Nconfs[i];j++)
	for(unsigned int k=0; k<raw_data[i][j].size();k++) {
	  stringstream s(raw_data[i][j][k]);
	  for(int col=0;col<Ncols;col++) {
	    double tmp;
	    s>>tmp;
	    convert_raw_data[col][i][j].push_back(tmp);
	  }
	}
  
  
    for(int col=0;col<Ncols;col++) {
      data.insert(pair<int,VVVfloat>(col, convert_raw_data[col]));
      for(int iens=0; iens< this->size; iens++) Transpose_VV(data[col][iens]);
    }

}





  void file_t::Read() {

    ifstream Read(path);

    if(!Read.is_open()) crash("In File_t Read() unable to open file: "+path);
    string ReadLine="";
    int it=0;
    int ncols_old;
    int Ncols;
    do {
      do {getline(Read, ReadLine);} while(ReadLine=="" && !Read.eof());
      int ncols=0;
      stringstream s(ReadLine);
      double tmp=0;
      while(s>>tmp) ncols++; 
      if(ncols != ncols_old && it!=0 && !Read.eof() && ReadLine.find(SC) == string::npos) {
	cout<<"Number of columns in file: "+path+" is not constant over configs"<<endl;
	cout<<"ncols_old: "<<ncols_old<<" ncols: "<<ncols<<endl;
	crash("Exiting...");
      }
      if(it==0) data.resize(ncols);
      if(!Read.eof() && ReadLine.find(SC) == string::npos) {
	Ncols=ncols;
	ncols_old=ncols;
	it++;
	nrows++;
      }
      if(!Read.eof() && ReadLine.find(SC) == string::npos) {
	stringstream ss(ReadLine);
	double tmp2;
	int c=0;
	while(ss>>tmp2) {data[c].push_back(tmp2);c++;} 

      }
    } while( !Read.eof());

    return;
  }




Vfloat Read_From_File(string path, int icol, int ncols, int rows_to_skip) {

  ifstream Read(path);
  Vfloat res;

 

  if(!boost::filesystem::is_regular_file(path)) crash("In Read_From_File cannot open file name: "+path);

  string fake_string;
  for(int iskip=0;iskip<rows_to_skip;iskip++) getline(Read, fake_string); 
  
  do {

    Vfloat Col_read;

    if(!Read.eof()) { 
    for(int i=0;i<ncols;i++) { double a; Read>>a; Col_read.push_back(a);}
    }

    if(!Read.eof())  res.push_back( Col_read[icol]);

  } while(!Read.eof());


  // printV(res,path,1);
  return res;

}
