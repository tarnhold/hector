/*! \file   MLEBase.cpp
 *  \author Machiel Bos
 *
 * The class that prepares the arrays before MLE can be computed using
 * AmmarGrag or FullCov.
 *
 * \date 13/ 1/2012   Coimbra Library
 * \date  7/10/2012   Santa Clara
 * \date  9/ 6/2015   Santa Clara
 */
//==============================================================================
  #include "MLEBase.h"
  #include <iostream>
  #include <ostream>
  #include <fstream>
  #include <iomanip>
  #include <cstring>
  #include <cstdlib>
  #include <cmath>


//==============================================================================
// Subroutines
//==============================================================================


//---!!---------------------------------------
  MLEBase::MLEBase(void) : NaN(sqrt(-1.0)),
                           tpi(8.0*atan(1.0))
//---!!---------------------------------------
  {
    Control        &control = Control::getInstance();
    Observations   &observations=Observations::getInstance();
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    NoiseModel     &noisemodel=NoiseModel::getInstance(); 

    (void) control;

    using namespace std;

    //--- The matrices H and F and vector x are not really needed in this
    //    class. However, they are needed in AmmarGrag.cpp which is reading
    //    the matrices from here! Think of MLEBase as BaseClass. 

    //---- Just be elegant
    x = H = NULL;

    observations.get_values(m,&t,&x);
    designmatrix->get_H(n,&H);
    
    Nparam = noisemodel.get_Nparam();
    Ngaps  = observations.number_of_gaps();
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
//-----------------------
  MLEBase::~MLEBase(void)
//-----------------------
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
//-------------------------------------
  void MLEBase::show_leastsquares(void)
//-------------------------------------
  {
    int            i;
    double         *error;
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    Observations   &observations=Observations::getInstance();

    using namespace std; 
    //---- Already show value of sigma_eta
    cout << "STD of the driving noise: " << sigma_eta << endl;
 
    error = new double[n];
    for (i=0;i<n;i++) error[i] = sigma_eta*sqrt(C_theta[i+n*i]);
    designmatrix->show_results(theta,error);
    designmatrix->compute_xhat(theta);
    observations.save_mom(true);

    delete[] error;
  }



/*! Computed the likelihood function. If quick==true then I assume noise
 *  parameters have not changed and that therefore there is no need to
 *  prepare the covariance (compute Cholesky or l1 & l2). One can directly
 *  compute Least-Squares.
 */
//---------------------------------------------------
  double MLEBase::compute(double *param_, bool quick)
//---------------------------------------------------
  {
    int           N;
    double        l,penalty,*param=NULL;
    NoiseModel    &noisemodel=NoiseModel::getInstance();

    using namespace std;
    if (quick==true) {
      compute_LeastSquares(param_);
      penalty = 0.0;
    } else {
      //--- Copy noise parameters to auxiliary array 'param'
      if (Nparam>0) {
        param = new double[Nparam];
        cblas_dcopy(Nparam,param_,1,param,1);
      }

      //--- Penalty function must be called first to make sure the noise
      //    parameters are brought back to valid values before computing a 
      //    covariance matrix.
      penalty = noisemodel.compute_penalty(param);

      //--- Compute WLS (virtual subroutine which is defined in the derived
      //    classes). This automatically sets the values of 'ln_det_C' and
      //    'sigma_eta'.
      prepare_covariance(param);
      compute_LeastSquares(param);
    }

    N = m - Ngaps;
    l = 0.5*(N*log(tpi) + penalty + ln_det_C + 2.0*N*log(sigma_eta) + N);

    if (param!=NULL) delete[] param;

    return l;
  }



/*! Compute log-likelihood and Information Criteria
 */
//----------------------------------------------
  void MLEBase::compute_L_and_ICs(double *param)
//----------------------------------------------
  {
    using namespace std;
    int             k,N;

    //--- How many parameters are estimated? 
    k = Nparam + n + 1;

    //--- Compute log(likelihood) at minimum
    N     = m - Ngaps;
    ln_L  = -compute(param,true);
    AIC   = 2.0*k - 2.0*ln_L;
    BIC   = k*log(N) - 2.0*ln_L;
  }



/*! For debugging it's convenient to have a simple way of showing a
 *  matrix on the screen (and to a file).
 */
//---------------------------------------------------------------------
  void MLEBase::show_matrix(const char name[], double *A, int m, int n)
//---------------------------------------------------------------------
  {
    using namespace std;
    fstream fp;
    char    filename[80];
    int     i,j;

    sprintf(filename,"%s.dat",name);
    fp.open(filename,ios::out);
    cout << name << "=" << endl;
    for (i=0;i<m;i++) {
      for (j=0;j<n;j++) printf("%7.3lf  ",A[i+m*j]);
      cout << endl;
      for (j=0;j<n;j++) fp << setprecision(5) << setw(10) << A[i+m*j] << "  ";
      fp << endl;
    }
    fp.close();
  }
