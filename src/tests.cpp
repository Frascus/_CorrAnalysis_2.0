#include "../include/tests.h"


void test() {




  int N=500;
  double texp=10;
  int R2=300;
  int Njacks=(int)(1.0*(N-R2)/(max(1.0,texp)));

  vector<int> counts;
  for(int i=0;i<N;i++) counts.push_back(i);

  Eigen::VectorXd params_val(N);
  for(int i=0;i<N;i++) params_val(i)=1;


  cout<<"Using Njacks: "<<Njacks<<endl;
 



  //Read covariance matrix
  Eigen::MatrixXd Cov(N, N);
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++) Cov(i,j) = exp(-abs(i-j)/texp);

  GaussianMersenne G(230765);
  Vfloat R=Covariate(Cov, params_val, G);
  


  

  //define std and improved estimator

  auto std_est= [&](const Vfloat &V) -> distr_t {

    int Nn= V.size();
    distr_t JACK(1);
    double bs = ((double)Nn)/((double)Njacks); //fractional block_size
 
    //get total sum
    double total_sum=0.0;
    for(int iconf=0;iconf< Nn;iconf++) total_sum += V[iconf];
    
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
        data_to_cluster -= V[iConf]*rect_size;



	//update position
	binPos = rpoint;


      } while (bin_end - binPos > 1e-10);
      

      data_to_cluster /= (double)(Nn - bs);

      JACK.distr.push_back(data_to_cluster);
    }
    
    
    return JACK;
  };


  auto impr_est= [&](const Vfloat &V, const vector<int> &dist) -> distr_t {

    distr_t JACK(1);
    int Nn=V.size();

    
    Vfloat Vnew;
    for(int i=0;i<Nn;i++) {

      Vnew.push_back(V[i]);
      
      for(int j=1;j<dist[i];j++) {
	if(dist[i]%2 != 0) { //d is odd
	  if( j<=(dist[i]-1)/2) Vnew.push_back(V[i]);
	  else Vnew.push_back( V[i+1]   );
	}
	else { //d is even
	  if(j<(dist[i]-1.0)/2.0) Vnew.push_back(V[i]);
	  else if (j==dist[i]/2) Vnew.push_back( 0.5*(V[i]+V[i+1]));
	  else Vnew.push_back(V[i+1]);
	}
      }
    }

    Nn=Vnew.size();
 
     //get total sum
    double total_sum=0.0;
    for(int iconf=0;iconf< Nn;iconf++) total_sum += Vnew[iconf];

    double bs = ((double)Nn)/((double)Njacks); //fractional block_size
    
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
        data_to_cluster -= Vnew[iConf]*rect_size;



	//update position
	binPos = rpoint;


      } while (bin_end - binPos > 1e-10);
      

      data_to_cluster /= (double)(Nn - bs);

      JACK.distr.push_back(data_to_cluster);
    }
    
    
    return JACK;
  };
  

  distr_t EST_ORIGINAL= std_est(R);
  cout<<"ORIGINAL: "<<EST_ORIGINAL.ave()<<" +- "<<EST_ORIGINAL.err()<<endl;


  //exclude R1 points from R  at random

  int Nreplica=1;
  RandomMersenne RM(54363463,R2);
  
  for(int irep=0;irep<Nreplica;irep++) {

    random_shuffle(counts.begin()+1,counts.end()-1);
    
    vector<int> R_to_erase=counts; R_to_erase.erase(R_to_erase.begin());  R_to_erase.resize(R2); sort(R_to_erase.begin(), R_to_erase.end());

    vector<int> dist(N-R2,1);
    for(int r=0;r<R2;r++) {
	dist[R_to_erase[r]-1-r]++;
    }
    //sort(counts.begin(), counts.end());
 
    Vfloat R_new= R;
    for(int r=0;r<R2;r++) R_new.erase( R_new.begin() + R_to_erase[R2-r-1]);

    distr_t EST_STD= std_est(R_new);
    distr_t EST_IMPR= impr_est(R_new, dist);

    cout<<"EST_STD: "<<EST_STD.ave()<<" +- "<<EST_STD.err()<<endl;
    cout<<"EST_IMPR: "<<EST_IMPR.ave()<<" +- "<<EST_IMPR.err()<<endl;
    //printV(R,"R", 0);
    //printV(R_to_erase, "R_to_erase", 0);
    //printV(R_new, "R_new",0);

  }
  
  
  


    
  return;
}
