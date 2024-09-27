#ifndef __highPrec__
#define __highPrec__

#include <mpfr.h>
#include "numerics.h"

using namespace std;

/// Undef this to disable actual usage of the high precision
//#define FAKE_HP

/// Structure to represent arbitray precision real number
struct PrecFloat
{
#ifndef FAKE_HP
  static inline int& _defaultPrecision()
  {
    static int _val;
    
    return _val;
  }
#endif
  
  /// Sets the default precision
  static void setDefaultPrecision(const int& n)
  {
#ifndef FAKE_HP
    _defaultPrecision()=n;
#endif
  }
  
  /// Gets the default precision
  static int getDefaultPrecision()
  {
    return
#ifndef FAKE_HP
      _defaultPrecision()
#else
      53
#endif
      ;
  }
  
  /// Gets the number of decimal digits that can be printed
  static int getNDigits()
  {
    return
#ifndef FAKE_HP
      getDefaultPrecision()*log10(2.0)
#else
      16
#endif
      ;
  }

  /// Returns the maximum exponent emax: ( max num is (1-e)x 2^emax )
  static double getEmax_max()
  {
    return
#ifndef FAKE_HP
    mpfr_get_emax_max()
#else
    0
#endif
      ;
  }

  /// Returns the current maximum exponent
  static double getEmax()
  {
    return
#ifndef FAKE_HP
    mpfr_get_emax()
#else
    0
#endif
      ;
  }
  
  /// Returns the current smaller number
  static PrecFloat getEpsilon()
  {
    return
      pow((PrecFloat)2,-getDefaultPrecision());
  }
  
  /// Storage
#ifndef FAKE_HP
  mpfr_t data{};
#else
  double data;
#endif
  
  inline friend std::ostream& operator<<(std::ostream& os,const PrecFloat& f)
  {
#ifndef FAKE_HP
    const bool fixed=
      os.flags()&std::ios_base::fixed;
       
    constexpr int lenFormat=10;
    char format[lenFormat];
    snprintf(format,lenFormat,"%%" "." "%td" "R%c",os.precision(),(fixed?'f':'g'));
    
    constexpr int lenOut=1024;
    char out[lenOut];
    mpfr_snprintf(out,lenOut,format,f.data);
    
    os<<out;
#else
    os<<f.data;
#endif
    
    return os;
  }


  
  /// Returns the internal data
  double get() const
  {
#ifdef FAKE_HP
    return data;
#else
    return
      mpfr_get_d(data,MPFR_RNDD);
#endif
  }
  
  /// Assignment
  PrecFloat& operator=(const double& in)
  {
#ifdef FAKE_HP
    data=in;
#else
    mpfr_set_d(data,in,MPFR_RNDD);
#endif
    
    return
      *this;
  }
  
  /// Assign from another number
  PrecFloat& operator=(const PrecFloat& oth)
  {
#ifdef FAKE_HP
    data=oth.data;
#else
    mpfr_set(data,oth.data,MPFR_RNDD);
#endif
    
    return
      *this;
  }
  
  /// Default initialization
  PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_init2(data,_defaultPrecision());
#endif
  }
  
  /// Copy constructor
  PrecFloat(const PrecFloat& oth) :
    PrecFloat()
  {
#ifdef FAKE_HP
    data=oth.data;
#else
    mpfr_set(data,oth.data,MPFR_RNDD);
#endif
  }
  
  /////////////////////////////////////////////////////////////////
  
#ifdef FAKE_HP
#define PROVIDE_CONVERSION_FROM(TYPE,MPFR_TAG)		\
  PrecFloat(const TYPE& in)				\
  {							\
    data=in;						\
  }
#else
#define PROVIDE_CONVERSION_FROM(TYPE,MPFR_TAG)		\
  PrecFloat(const TYPE& in) :				\
    PrecFloat()						\
  {							\
    mpfr_set_ ## MPFR_TAG(data,in,MPFR_RNDD);		\
  }
#endif
  
  PROVIDE_CONVERSION_FROM(double,d);
  PROVIDE_CONVERSION_FROM(int,si);
  PROVIDE_CONVERSION_FROM(unsigned long int,ui);
  
#undef PROVIDE_CONVERSION_FROM
  
  /// Destructor
  ~PrecFloat()
  {
#ifdef FAKE_HP
#else
    mpfr_clear(data);
#endif
  }
  
  /////////////////////////////////////////////////////////////////
  
#ifdef FAKE_HP
#define BINARY_HELPER(NAME,MPFR_NAME)		\
  out.data=NAME(in1.data,in2.data)
#else
#define BINARY_HELPER(NAME,MPFR_NAME)			\
  MPFR_NAME(out.data,in1.data,in2.data,MPFR_RNDD)
#endif
  
#define PROVIDE_BINARY_FUNCTION(NAME,MPFR_NAME)		\
  							\
  inline friend PrecFloat NAME(const PrecFloat& in1,	\
			       const PrecFloat& in2)	\
  {							\
    PrecFloat out;					\
    							\
    BINARY_HELPER(NAME,MPFR_NAME);			\
    							\
    return						\
      out;						\
  }
  
  PROVIDE_BINARY_FUNCTION(pow,mpfr_pow)
  
#undef PROVIDE_BINARY_FUNCTION
#undef BINARY_HELPER
  
  /////////////////////////////////////////////////////////////////
  
#ifdef FAKE_HP
#define BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME)	\
  data NAME in.data
#else
#define BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME)	\
  MPFR_NAME(data,in.data)
#endif
  
#define PROVIDE_BINARY_COMPARISON_OPERATOR(NAME,MPFR_NAME)	\
  								\
  inline bool operator NAME(const PrecFloat& in) const	\
  {								\
    return							\
      BINARY_COMPARISON_OPERATOR_HELPER(NAME,MPFR_NAME);	\
  }

  PROVIDE_BINARY_COMPARISON_OPERATOR(<,mpfr_less_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(<=,mpfr_lessequal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(>,mpfr_greater_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(>=,mpfr_greaterequal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(==,mpfr_equal_p);
  PROVIDE_BINARY_COMPARISON_OPERATOR(!=,!mpfr_equal_p);
  
#undef BINARY_COMPARISON_OPERATOR_HELPER
#undef BINARY_COMPARISON_OPERATOR
  
  /////////////////////////////////////////////////////////////////
  
  // Providing the binary operator as nonmember function allows to
  // take into account automatically also the cases in which the first
  // operand is not a PrecFloat
  
#ifdef FAKE_HP
#define BINARY_OPERATOR_HELPER(NAME,MPFR_NAME)\
  out.data=in1.data NAME in2.data
#else
#define BINARY_OPERATOR_HELPER(NAME,MPFR_NAME)\
  MPFR_NAME(out.data,in1.data,in2.data,MPFR_RNDD)
#endif
  
#define PROVIDE_SELF_BINARY_OPERATOR(NAME)			\
  								\
  inline PrecFloat& operator NAME ## =(const PrecFloat& in)	\
    {								\
      return							\
	(*this)=(*this)NAME in;					\
    }
  
#define PROVIDE_BINARY_OPERATOR(NAME,MPFR_NAME)			\
  								\
  friend inline PrecFloat operator NAME(const PrecFloat& in1,	\
					const PrecFloat& in2)	\
  {								\
    PrecFloat out;						\
  								\
    BINARY_OPERATOR_HELPER(NAME,MPFR_NAME);			\
								\
    return out;							\
  }								\
  								\
  PROVIDE_SELF_BINARY_OPERATOR(NAME);				\
  
  PROVIDE_BINARY_OPERATOR(+,mpfr_add)
  PROVIDE_BINARY_OPERATOR(-,mpfr_sub)
  PROVIDE_BINARY_OPERATOR(*,mpfr_mul)
  PROVIDE_BINARY_OPERATOR(/,mpfr_div)
  
#undef BINARY_OPERATOR_HELPER
#undef PROVIDE_BINARY_OPERATOR
#undef PROVIDE_SELF_BINARY_OPERATOR
  
  /// Negation
  PrecFloat operator-() const
  {
    PrecFloat out;
    
#ifdef FAKE_HP
    out.data=-data;
#else
    mpfr_neg(out.data,this->data,MPFR_RNDD);
#endif
    
    return
      out;
  }
};

/////////////////////////////////////////////////////////////////

#ifdef FAKE_HP
#define UNARY_HELPER(NAME,MPFR_NAME)\
  out.data=NAME(in.data)
#else
#define UNARY_HELPER(NAME,MPFR_NAME)\
  MPFR_NAME(out.data,in.data,MPFR_RNDD)
#endif

#define PROVIDE_UNARY_FUNCTION(NAME,MPFR_NAME)	\
						\
  inline PrecFloat NAME(const PrecFloat& in)	\
{						\
  PrecFloat out;				\
						\
  UNARY_HELPER(NAME,MPFR_NAME);			\
						\
  return					\
    out;					\
}

PROVIDE_UNARY_FUNCTION(exp,mpfr_exp)
PROVIDE_UNARY_FUNCTION(abs,mpfr_abs)
PROVIDE_UNARY_FUNCTION(sqrt,mpfr_sqrt)
PROVIDE_UNARY_FUNCTION(asin,mpfr_asin)
PROVIDE_UNARY_FUNCTION(acos,mpfr_acos)
PROVIDE_UNARY_FUNCTION(atan,mpfr_atan)
PROVIDE_UNARY_FUNCTION(sin,mpfr_sin)
PROVIDE_UNARY_FUNCTION(cos,mpfr_cos)
PROVIDE_UNARY_FUNCTION(tan,mpfr_tan)
PROVIDE_UNARY_FUNCTION(sinh,mpfr_sinh)
PROVIDE_UNARY_FUNCTION(cosh,mpfr_cosh)
PROVIDE_UNARY_FUNCTION(tanh, mpfr_tanh)
PROVIDE_UNARY_FUNCTION(sech, mpfr_sech)
PROVIDE_UNARY_FUNCTION(erf,mpfr_erf)
PROVIDE_UNARY_FUNCTION(erfc,mpfr_erfc)
PROVIDE_UNARY_FUNCTION(log, mpfr_log)
PROVIDE_UNARY_FUNCTION(gamma, mpfr_gamma)


#undef PROVIDE_UNARY_FUNCTION
#undef UNARY_HELPER

#ifdef FAKE_HP
#define BINARY_HELPER(NAME, MPFR_NAME)		\
  out.data=NAME(in1.data, in2.data)
#else
#define BINARY_HELPER(NAME,  MPFR_NAME)			\
  MPFR_NAME(out.data, in1.data, in2.data, MPFR_RNDD)
#endif

#define PROVIDE_BINARY_FUNCTION(NAME, MPFR_NAME)			\
  									\
  inline PrecFloat NAME(const PrecFloat& in1, const PrecFloat& in2)	\
  {									\
    PrecFloat out;							\
    									\
    BINARY_HELPER(NAME, MPFR_NAME);					\
    									\
    return								\
      out;								\
  }

PROVIDE_BINARY_FUNCTION(min, mpfr_min)
PROVIDE_BINARY_FUNCTION(max, mpfr_max)
PROVIDE_BINARY_FUNCTION(gamma_inc, mpfr_gamma_inc)



#undef PROVIDE_BINARY_FUNCTION
#undef BINARY_HELPER


/// Returns the square
inline PrecFloat sqr(const PrecFloat& in)
{
  return in*in;
}

/////////////////////////////////////////////////////////////////

/// Precise definition of Pi
inline PrecFloat precPi()
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=M_PI;
#else
  mpfr_const_pi(out.data,MPFR_RNDD);
#endif
  
  return
    out;
}


/// Precise Euler constant

inline PrecFloat precEuler()
{
  PrecFloat out;
  
#ifdef FAKE_HP
  out.data=0.57721566490153286060651209008240243104215933593992;
#else
  mpfr_const_euler(out.data, MPFR_RNDD);
#endif
  
  return
    out;
}


/////////////////////////////////////////////////////////////////

//return factorial of no
inline PrecFloat precFact(unsigned long int n)
{
  PrecFloat out;

#ifdef FAKE_HP
  out.data= fact(n);
#else
  mpfr_fac_ui(out.data,n, MPFR_RNDD); 
#endif

  return
    out;
}

// Tell Eigen how to deal with PrecFloat numbers

namespace Eigen
{
  template<>
  struct NumTraits<PrecFloat> :
    GenericNumTraits<PrecFloat>
  {
    typedef PrecFloat Real;
    typedef PrecFloat NonInteger;
    typedef PrecFloat Nested;
    
    static inline Real epsilon()
    {
      return 0;
    }
    
    static inline Real dummy_precision()
    {
      return 0;
    }
    
    static inline int digits10()
    {
      return 0;
    }
    
    enum
      {
	IsInteger=0,
	IsSigned=1,
	IsComplex=0,
	RequireInitialization=1,
	ReadCost=6,
	AddCost=150,
	MulCost=100
      };
  };
}

using PrecVect=
  Eigen::Matrix<PrecFloat,Eigen::Dynamic,1>;

using PrecMatr=
  Eigen::Matrix<PrecFloat,Eigen::Dynamic,Eigen::Dynamic>;

/// Integrate the passed function f with the double exponential
/// transformation of eq.1.16 of
/// https://core.ac.uk/download/pdf/82080804.pdf
///
/// The function c takes into account the jacobian and the change of variable
///
/// The trapizio integration is carried out recursively, halving the
/// stepsize until the maximal relative precision is reached.
///
/// The value at zero is evaluated before everything and used as a
/// first estimate of the integral, I0. Then, at each iteration i, the
/// value I(i-1) is considered as \sum_i=0^n c(x*2*i*step)*2*step
/// Then, the estimate I(i) is obtained as
/// I(i) = I(i-1)/2 + \sum_i=0^n c(x*(2*i+1)*step)*step
///
/// n is taken to be the first point where the term is negligible
///
/// n might have to be adapted during the i-itertion, in case the
/// needed part of the full summation is re-evaluated
///
/// The iterative procedure is stopped when the result is stable
/// within attainable precisione
template <typename F>
PrecFloat integrateUpToInfinite(F&& f,const double& xMin=0.0,const int& verbose=0)
{

  int numDigits= PrecFloat::getNDigits();
  int original_precision = PrecFloat::getDefaultPrecision();
  auto old_cout_precision= cout.precision( numDigits);

  
  const PrecFloat piHalf=
    precPi()/2;
 

  if(verbose) cout<<"integrateUpToInfinite called with: xMin: "<<xMin<<endl;
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  /// We compute up to the highest precision possible, which needs to
  /// be adjusted in terms of the number of iterations (maybe it might
  /// be enough to increase with the square root?)

  PrecFloat InitialmaxAttainableStability= PrecFloat::getEpsilon()*10;

  PrecFloat sum_final;
  bool CHANGE_PRECISION=false;
  int i=0;
  
  do {

    CHANGE_PRECISION=false;  
    int COUNT_RESTART=0;
    if(verbose) cout<<"Precision set to : "<<(original_precision+original_precision*i/2)<<endl<<flush;
    PrecFloat::setDefaultPrecision( original_precision + original_precision*i/2);
    PrecFloat maxAttainableStability= InitialmaxAttainableStability;
  
  
   
    auto c=
      [&f,&xMin, &piHalf, &verbose](const PrecFloat& t)
      {
	const PrecFloat s=sinh(t);
	const PrecFloat x=exp(piHalf*s)+xMin;
	const PrecFloat fx= f(x);     
	const PrecFloat jac=piHalf*exp(piHalf*s)*cosh(t);
	const PrecFloat res=fx*jac; 

	/*
	if(verbose) {
	  cout<<"##########"<<endl;
	  cout<<"t: "<<t<<" s: "<<sinh(t)<<" x: "<<x<<endl;
	  cout<<"x: "<<x<<" f(x): "<<fx<<", jac(x): "<<jac<<endl;
	}
	*/
	
      
	return res;
      };

    PrecFloat step_old=1;
    bool RESTART=false;
    PrecFloat stability;
    
    do {
      int it=0;
      RESTART=false;  
      PrecFloat sum=c(0)*2*step_old;
      PrecFloat extreme=0;
      PrecFloat step=step_old;
      PrecFloat precSum;
     
      do
	{
	  precSum=sum;
	  
	  bool converged=
	    false;
	  
	  //bool converged_minus=false;
	  //bool converged_plus=false;
	  
	  sum/=2;
	  PrecFloat t=step;
	  
	  const PrecFloat doubleStep=
	    step*2;
	  
	  bool exitTheLoop=
	    false;

	  int nsteps=0;
	  
	  while(not exitTheLoop)
	    {
	      nsteps++;
	      const PrecFloat contr=c(t)+c(-t);
	      PrecFloat newSum= sum + contr*step;
	    
	      //converged_plus= ( abs(contr_plus*step/sum) < PrecFloat(0.1)*PrecFloat::getEpsilon());
	      //converged_minus= ( abs(contr_minus*step/sum) < PrecFloat(0.1)*PrecFloat::getEpsilon());
	      
	      converged= ( abs(contr*step/(sum)) < PrecFloat(0.1)*PrecFloat::getEpsilon());
	      //converged = (sum==newSum);
	      
	      if(verbose) {
		cout<<"RESTART: "<<COUNT_RESTART<<" t: "<<t<<" step: "<<step<<" contr: "<<contr<<" t>extreme: "<<(t>extreme)<<" extreme: "<<extreme<<" converged: "<<converged<<endl<<flush;
		cout<<"LOOPING: sum: "<<newSum<<" precSum: "<<sum<<endl<<flush;
		cout.precision( 2*numDigits);
		cout<<"LOOPING(x2 prec): sum: "<<newSum<<" precSum: "<<sum<<endl<<flush;
		cout.precision( numDigits);
	      }
	      
		  
	
	  
	      exitTheLoop=
		(converged and t>extreme);

	    	  
	      if(t>extreme)
		{
		  extreme=t;
		  t+=step;
		}
	      else
		t+=doubleStep;
	  
	      sum=newSum;
	    };
	
	  step/=2;
	

	  
	  stability=abs(sum/precSum-1);
	  
	  if(verbose) {
	    cout<<"LOOP EXITED: sum: "<<sum<<" precSum: "<<precSum<<", extreme: "<<extreme<<" step: "<<step*2<<endl<<flush;
	    cout<<"Stability: "<<stability<<" MaxAttainable: "<<maxAttainableStability<<endl<<flush;
	  }
	  
	  maxAttainableStability*=2;

	  if( stability <= maxAttainableStability) {sum_final=sum; RESTART=false; CHANGE_PRECISION=false;}
	  else {
	    if(COUNT_RESTART==0) {
	       if(step_old/step > PrecFloat(5e3)) {RESTART=true; step_old=step; CHANGE_PRECISION=false;}
	    }
	    else {
	      if( (COUNT_RESTART >=1) && (it >=5) ) CHANGE_PRECISION=true;
	    }
	  }

	  if(RESTART) COUNT_RESTART++;
	 
	  
	
	  it++;
	}
      while( (stability>maxAttainableStability) && (!RESTART)  && (!CHANGE_PRECISION) );

      if(verbose) {
	if(RESTART) cout<<"Restarting..."<<endl<<flush;
	else cout<<"Restart not needed!"<<endl<<flush; }
    }
    while(RESTART && !CHANGE_PRECISION);
    i++;

    if(CHANGE_PRECISION && verbose) cout<<"Changing precision"<<endl<<flush;
    if(!CHANGE_PRECISION && verbose) cout<<"Attained precision: "<<stability<<endl<<flush;
    
  } while( CHANGE_PRECISION );
  
  PrecFloat::setDefaultPrecision(original_precision);
  
  if(verbose) {
    cout<<"final sum: "<<sum_final<<endl<<flush;
    cout<<"integral computed"<<endl<<flush;
  }
  
  if(isnan(sum_final.get())) crash("In integrateUpToInfinity res is nan");

  cout.precision(old_cout_precision);

   
  return sum_final;
}




template <typename F>
PrecFloat integrateUpToXmax(F&& f,const double& xMin=0, const double& xMax=1,const int& verbose=0)
{

  int numDigits= PrecFloat::getNDigits();
  int original_precision = PrecFloat::getDefaultPrecision();
  auto old_cout_precision= cout.precision( numDigits);



  PrecFloat piHalf= precPi()/2;
  
  if(verbose) cout<<"integrateUpToXmax called with: xMin: "<<xMin<<endl;
  
  /// We compute up to the highest precision possible, which needs to                                                                                                                                    
  /// be adjusted in terms of the number of iterations (maybe it might                                                                                                                                   
  /// be enough to increase with the square root?)                                                                                                                                                       

  PrecFloat InitialmaxAttainableStability= PrecFloat::getEpsilon()*10;

  PrecFloat sum_final;
  bool CHANGE_PRECISION=false;
  int i=0;

  do {

    CHANGE_PRECISION=false;
    int COUNT_RESTART=0;
    if(verbose) cout<<"Precision set to : "<<(original_precision+original_precision*i/2)<<endl<<flush;
    PrecFloat::setDefaultPrecision( original_precision + original_precision*i/2);
    PrecFloat maxAttainableStability= InitialmaxAttainableStability;
    
    
   
    auto c=
      [&f, &xMin, &xMax, &piHalf, &verbose](const PrecFloat& t)
      {
	//tanh-sinh to map [xMin, xMax] into [-infinite, +infinite]
	PrecFloat xMin_prec=xMin; PrecFloat xMax_prec=xMax;
      
	const PrecFloat s=sinh(t);
	const PrecFloat g=tanh(piHalf*s);
	const PrecFloat Cg= sech(piHalf*s);
	const PrecFloat x=g*(xMax_prec-xMin_prec)/2+(xMax_prec +xMin_prec)/2;
	const PrecFloat jac=((xMax_prec-xMin_prec)/2)*piHalf*Cg*Cg*cosh(t);
	const PrecFloat fx= f(x);

	if(verbose) {
	  cout<<"##########"<<endl;
	  cout<<"t: "<<t<<" s: "<<sinh(t)<<" x: "<<x<<endl;
	  cout<<"x: "<<x<<" f(x): "<<fx<<", jac(x): "<<jac<<endl;
	}

	
	return fx*jac;

      };
  
    PrecFloat step_old=1;
    PrecFloat stability;
    bool RESTART=false;

    do {

      int it=0;
      RESTART=false;  
      PrecFloat sum=c(0)*2*step_old;
      PrecFloat extreme=0;
      PrecFloat step=step_old;
      PrecFloat precSum;
   
 
      do
	{
	  precSum=sum;
      
	  bool converged=
	    false;
      
	  sum/=2;
	  PrecFloat t=step;
      
	  const PrecFloat doubleStep=
	    step*2;
	
	  bool exitTheLoop=
	    false;
      
	  while(not exitTheLoop)
	    {
	      const PrecFloat contr=
		c(t)+c(-t);
	  
	      PrecFloat newSum= sum + contr*step;
		
	      converged= ( abs(contr*step/sum) < PrecFloat(0.1)*PrecFloat::getEpsilon());

	      
	      if(verbose) {
		cout<<"RESTART: "<<COUNT_RESTART<<" t: "<<t<<" step: "<<step<<" contr: "<<contr<<" t>extreme: "<<(t>extreme)<<" extreme: "<<extreme<<" converged: "<<converged<<endl<<flush;
		cout<<"LOOPING: sum: "<<newSum<<" precSum: "<<sum<<endl<<flush;
		cout.precision( 2*(PrecFloat::getNDigits()));
		cout<<"LOOPING(x2 prec): sum: "<<newSum<<" precSum: "<<sum<<endl<<flush;
		cout.precision( PrecFloat::getNDigits());
	      }
	      
	  
	      exitTheLoop=
		(converged and t>extreme);
	  
	      if(t>extreme)
		{
		  extreme=t;
		  t+=step;
		}
	      else
		t+=doubleStep;
	      
	      sum=newSum;
	    };
	  
	  step/=2;
	
	
      
	  stability=abs(sum/precSum-1);

	  if(verbose) {
	    cout<<"LOOP EXITED: sum: "<<sum<<" precSum: "<<precSum<<", extreme: "<<extreme<<" step: "<<step*2<<endl<<flush;
	    cout<<"Stability: "<<stability<<" MaxAttainable: "<<maxAttainableStability<<endl<<flush;
	  }
	  
	  maxAttainableStability*=2;
	  
	  if( stability <= maxAttainableStability) {sum_final=sum; RESTART=false; CHANGE_PRECISION=false;}
	  else {
	    if(COUNT_RESTART==0) {
	      if(step_old/step > PrecFloat(5e3)) {RESTART=true; step_old=step; CHANGE_PRECISION=false;}
	    }
	    else {
	      if( (COUNT_RESTART >=1) && (it >=5) ) CHANGE_PRECISION=true;
	    }
	  }

	  if(RESTART) COUNT_RESTART++;

	  it++;

		   
	}
      while( (stability>maxAttainableStability) && (!RESTART) && (!CHANGE_PRECISION)   );

      if(verbose) {
	if(RESTART) cout<<"Restarting..."<<endl<<flush;
	else cout<<"Restart not needed!"<<endl<<flush; }
    }
    while(RESTART && !CHANGE_PRECISION);

    i++;

    if(CHANGE_PRECISION && verbose) cout<<"Changing precision"<<endl<<flush;
    if(!CHANGE_PRECISION && verbose) cout<<"Attained precision: "<<stability<<endl<<flush;

  } while( CHANGE_PRECISION );

  PrecFloat::setDefaultPrecision(original_precision);


  if(verbose) {
    cout<<"final sum: "<<sum_final<<endl<<flush;
    cout<<"integral computed"<<endl<<flush;
  }
  
  if(isnan(sum_final.get())) crash("In integrateUpToXmax res is nan");
  
  cout.precision(old_cout_precision);


  
  
  return sum_final;
}







PrecFloat ExpEiComplexSum(PrecFloat MOD, PrecFloat PH, PrecFloat s,  bool MODE); 
PrecFloat Erfi( PrecFloat x);
PrecFloat DawsonF( PrecFloat x);


#ifdef MAIN
int PrecFloat::_defaultPrecision;
#endif


#endif
