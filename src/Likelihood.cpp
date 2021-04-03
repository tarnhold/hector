/*! \file   Likelihood.cpp
 *  \author Machiel Bos
 *
 * The core class that computes the likelihood function. The exact method 
 * for computing the likelihood and Least-Squares bit is define in the derived
 * classes LevinsonGap, LevinsonNoGap, AmmarGrag and FullCov classes.
 *
 * \date 13/ 1/2012   Coimbra Library
 * \date  7/10/2012   Santa Clara
 */
//==============================================================================
  #include "Likelihood.h"
  #include <iostream>
  #include <ostream>
  #include <fstream>
  #include <iomanip>
  #include <cstring>
  #include <cstdlib>
  #include <cmath>

  //--- To avoid circular include problems, I put the following includes here
  //    The explanation is that Likelihood.h defines __LIKELIHOOD and then it
  //    includes "LevinsonGap.h" which again tries to include "Likelihood.h"
  //    because __LIKELIHOOD is now defined, the class definition is not
  //    read and the LevinsonGap class cannot determine well its base class,
  //    causing problems.
  #include "LevinsonNoGap.h"
  #include "LevinsonGap.h"
  #include "AmmarGrag.h"
  #include "FullCov.h"

//==============================================================================
// Subroutines
//==============================================================================

//++++++ Singleton stuff +++++++++++++++++
bool Likelihood::instanceFlag = false;
Likelihood* Likelihood::singleton = NULL;
//++++++++++++++++++++++++++++++++++++++++


//--!!-------------------------------------
  Likelihood* Likelihood::getInstance(void)
//--!!-------------------------------------
  {
    int           Ngaps;
    Control       *control = Control::getInstance();
    Observations  *observations = Observations::getInstance();
    char          methodname[80];

    using namespace std;
    Ngaps = observations->number_of_gaps();
    if (instanceFlag==false) {
      instanceFlag=true;
      control->get_string("LikelihoodMethod",methodname);
      cout << endl << "----------------" << endl << "  " << methodname
           << endl << "----------------" << endl;
      if (strcmp(methodname,"Levinson")==0) {
        if (Ngaps==0)  singleton = new LevinsonNoGap();
        else           singleton = new LevinsonGap();
      } else if (strcmp(methodname,"AmmarGrag")==0) {
        singleton = new AmmarGrag();
      } else if (strcmp(methodname,"FullCov")==0) {
        singleton = new FullCov();
      } else {
        cerr << "Unkown Likelihood Method: " << methodname << endl;
        exit(EXIT_FAILURE);
      }
    }
    return singleton;
  }



//---!!--------------------------------------------
  Likelihood::Likelihood(void) : tpi(8.0*atan(1.0))
//---!!--------------------------------------------
  {
    Observations   *observations=Observations::getInstance();
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    NoiseModel     *noisemodel=NoiseModel::getInstance();

    //---- Just be elegant
    x = H = NULL;

    observations->get_values(m,&t,&x);
    designmatrix->get_H(n,&H);
    
    Nparam = noisemodel->get_Nparam();
    Ngaps  = observations->number_of_gaps();
    //--- Create F-matrix if Ngaps>0.
    if (Ngaps>0) {
      designmatrix->get_F(&F);
    } else {
      F = NULL;
    }

    //--- For Least-Squares we need theta and C_theta. Note that C_theta and
    //    theta are set/filled by the derived classes. I have put them in
    //    this base clase because their values need to be shown in the end
    //    to the user.
    theta   = new double[n];
    C_theta = new double[n*n];
    memset(theta,0,(n)*sizeof(double));     // put theta to zero
    memset(C_theta,0,(n*n)*sizeof(double)); // put C_theta to zero
  }



/*! Free up Memory 
 */
//-----------------------------
  Likelihood::~Likelihood(void)
//-----------------------------
  {
    if (t!=NULL)       delete[]  t;
    if (x!=NULL)       delete[]  x;
    if (H!=NULL)       delete[]  H;
    if (F!=NULL)       delete[]  F;
    if (theta!=NULL)   delete[]  theta;
    if (C_theta!=NULL) delete[]  C_theta;
  }



/*!
 */
//----------------------------------------
  void Likelihood::show_leastsquares(void)
//----------------------------------------
  {
    int            i;
    double         *error;
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    Observations   *observations=Observations::getInstance();

    using namespace std; 
    //---- Already show value of sigma_eta
    cout << "STD of the innovation noise: " << sigma_eta << endl;
 
    error = new double[n];
    for (i=0;i<n;i++) error[i] = sigma_eta*sqrt(C_theta[i+n*i]);
    designmatrix->show_results(theta,error);
    designmatrix->compute_xhat(theta);
    observations->save_mom(true);

    delete[] error;
  }



/*! Computed the likelihood function. 
 */
//------------------------------------------
  double Likelihood::compute(double *param_)
//------------------------------------------
  {
    int           N;
    double        lndeterminant,l,penalty,*param=NULL;
    NoiseModel    *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- Copy noise parameters to auxiliary array 'param'
    if (Nparam>0) param = new double[Nparam];
    cblas_dcopy(Nparam,param_,1,param,1);

    //--- Penalty function must be called first to make sure the noise
    //    parameters are brought back to valid values before computing a 
    //    covariance matrix.
    penalty = noisemodel->compute_penalty(param);

    //--- Compute WLS (virtual subroutine which is defined in the derived
    //    classes).
    compute_LeastSquares(param,lndeterminant,sigma_eta);

    N = m - Ngaps;
    l = 0.5*(N*log(tpi) + penalty + lndeterminant + 2.0*N*log(sigma_eta) + N);

    if (param!=NULL) delete[] param;

    return l;
  }



/*! For debugging it's convenient to have a simple way of showing a
 *  matrix on the screen (and to a file).
 */
//------------------------------------------------------------------------
  void Likelihood::show_matrix(const char name[], double *A, int m, int n)
//------------------------------------------------------------------------
  {
    using namespace std;
    fstream fp;
    char    filename[80];
    int     i,j;

    sprintf(filename,"%s.dat",name);
    fp.open(filename,ios::out);
    cout << name << "=" << endl;
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) cout << A[i+m*j] << "  ";
      cout << endl;
      for (j=0;j<n;j++) fp << setprecision(4) << setw(9) << A[i+m*j] << "  ";
      fp << endl;
    }
    fp.close();
  }



/*! For debugging it's convenient to have a simple way of showing a
 *  Fourier transformed matrix on the screen (and to a file).
 */
//------------------------------------------------------------------------
  void Likelihood::show_Fmatrix(const char name[], fftw_complex *A, 
							     int m, int n)
//------------------------------------------------------------------------
  {
    using namespace std;
    int     i,j;
    fstream fp;
    char    filename[80];
    
    sprintf(filename,"%s.dat",name);
    fp.open(filename,ios::out);
    cout << name << "=" << endl;
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) cout << A[i+m*j][0] << " ";
      cout << endl; 
      for (j=0;j<n;j++) cout << A[i+m*j][1] << " ";
      cout << endl; 
      for (j=0;j<n;j++) fp << setprecision(4) << setw(9) << A[i+m*j][0] << "  ";
      fp << endl; 
      for (j=0;j<n;j++) fp << setprecision(4) << setw(9) << A[i+m*j][1] << "  ";
      fp << endl; 
    }
    fp.close();
  }
