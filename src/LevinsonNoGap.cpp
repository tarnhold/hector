/*! \file   LevinsonNoGap.cpp
 *  \author Machiel Bos
 *
 * A class that computes the likelihood function, including the Least-Squares
 * part, using the Levinson method to perform a Cholesky decomposition of the
 * Covariance matric. This is my classic approach, assuming that there are
 * no gaps.
 *
 * \date 25/1/2012  Santa Clara
 */
//==============================================================================
  #include "LevinsonNoGap.h"
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstdio>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!-----------------------------
  LevinsonNoGap::LevinsonNoGap(void)
//---!!-----------------------------
  {
    char          NoTrans='N';
    int           one=1,info;
    NoiseModel    *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- Sanity check
    if (Ngaps!=0) {
      cerr << "LevinsonNoGaps: Ngaps=" << Ngaps << "!!!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- get number of noise model parameters to estimate
    Nparam = noisemodel->get_Nparam();

    //--- From the base class Likelihood, the values of n & m are known.
    try {
      gamma_x = new double[m];
      r       = new double[m];
      y       = new double[m];
      A       = new double[m*n];
      dummyH  = new double[m*n];
      dummyx  = new double[m];
      U       = new double[2*m];
      V       = new double[2*m];
      dummy   = new double[m*n];
    }
    catch (bad_alloc) {
      cerr << "Levinson: Need more memory for variables!" << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Free up memory
 */
//---!!------------------------------
  LevinsonNoGap::~LevinsonNoGap(void)
//---!!------------------------------
  {
    if (gamma_x!=NULL)   delete[]  gamma_x;
    if (r!=NULL)         delete[]  r;
    if (y!=NULL)         delete[]  y;
    if (A!=NULL)         delete[]  A;
    if (dummy!=NULL)     delete[]  dummy;
    if (dummyH!=NULL)    delete[]  dummyH;
    if (dummyx!=NULL)    delete[]  dummyx;
    if (U!=NULL)         delete[]  U;
    if (V!=NULL)         delete[]  V;
  }



/*! Mactrick
 *
 * \param[in]  gamma_x       : autocovariance matrix of noise
 * \param[out] lndeterminant : logarithm of determinant of C
 * \param[out] A             : U*H = A where U'*U=Cinv
 * \param[out] y             : U*x = y
 */
//---------------------------------------------------------------------------
  void LevinsonNoGap::mactrick(double *gamma_x,
                    		double& lndeterminant_, double *A, double *y)
//---------------------------------------------------------------------------
  {
    double   sin_theta,cos_theta;
    int      k,k_old=0,k_new=1;

    using namespace std;
    //--- define the generators u and v
    cblas_dcopy(m,gamma_x,1,U,1);
    cblas_dscal(m,1.0/sqrt(gamma_x[0]),U,1);
    cblas_dcopy(m,U,1,V,1);
    V[0] = 0.0;

    //--- First determinant value (1st element on diagonal)
    lndeterminant_ = log(U[0]);

    //--- First solution vector values
    cblas_dcopy(n,H,m,A,m);
    cblas_dscal(n,1.0/U[0],A,m);
    y[0] = x[0]/U[0];
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        m-1,n,1,1.0,U+1,m,A,m,0.0,dummyH+1,m);
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        m-1,1,1,1.0,U+1,m,y,m,0.0,dummyx+1,m);
    for (k=1;k<m;k++) {
      if (U[(k-1) + m*k_old]<1.0e-4) {
        cout << "mactrick: small U[k-1 + m*k_old]=" << U[k-1 + m*k_old] << endl;
        exit(EXIT_FAILURE);
      }
      sin_theta = V[k + m*k_old]/U[(k-1) + m*k_old];
      cos_theta = sqrt(1.0-sin_theta*sin_theta);
      if (cos_theta<1.0e-4) {
        cout << "mactrick: cos_theta=" << cos_theta << endl;
        exit(EXIT_FAILURE);
      }

      cblas_dcopy(m-k,&V[k + m*k_old],1,&V[k + m*k_new],1);
      cblas_daxpy(m-k,-sin_theta,&U[k-1 + m*k_old],1,&V[k + m*k_new],1);
      cblas_dscal(m-k,1.0/cos_theta,&V[k + m*k_new],1);

      cblas_dcopy(m-k,&U[(k-1) + m*k_old],1,&U[k + m*k_new],1);
      cblas_dscal(m-k,cos_theta,&U[k + m*k_new],1);
      cblas_daxpy(m-k,-sin_theta,&V[k + m*k_new],1,&U[k + m*k_new],1);

      //--- Update determinant
      lndeterminant_ += log(U[k + m*k_new]);

      //--- Extend back-substitution
      cblas_dcopy(n,H+k,m,A+k,m);
      cblas_daxpy(n,-1.0,dummyH+k,m,A+k,m);
      cblas_dscal(n,1.0/U[k + m*k_new],A+k,m);
      y[k] = (x[k]-dummyx[k])/U[k + m*k_new];
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                m-k-1,n,1,1.0,&U[k+1 + m*k_new],m,A+k,m,0.0,dummy+(k+1),m);
      cblas_daxpy(n*m,1.0,dummy,1,dummyH,1);
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                m-k-1,1,1,1.0,&U[k+1 + m*k_new],m,y+k,m,0.0,dummy+(k+1),m);
      cblas_daxpy(m,1.0,dummy,1,dummyx,1);
      k_old = k_new;
      k_new = 1 - k_new;
    }

    lndeterminant_ *= 2.0;  // remember C=U'*U !!
  }



/*! Compute the Least-Squares bit to estimate theta and sigma_eta. At
 *  the same time determine the logarithm of the determinant of C. Note
 *  that sigma_eta is estimated from the whitened residuals, not using MLE.
 *
 * \param[in]  param          : array with noise parameters
 * \param[out] lndeterminant  : logarithm of the likelihood
 * \param[out] sigma_eta      : STD of whitened residuals
 */
//-------------------------------------------------------------------------
  void LevinsonNoGap::compute_LeastSquares(double *param, double& 
					lndeterminant_, double& sigma_eta_)
//-------------------------------------------------------------------------
  {
    int         info,i,j;
    double      product;
    NoiseModel  *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- I assume the noise parameters, and thus the covariance matrix, have 
    //    been changed. As a result, matrix A and y need again to be computed.
    cblas_dcopy(m*n,H,1,A,1);
    cblas_dcopy(m,x,1,y,1);

    //--- Create the 1st column of the Covariance matrix. Since this is
    //    a Toeplitz matrix, this column vector is sufficient.
    //    Note that this is for unit variance value of sigma_eta.
    noisemodel->get_covariance(param,m,gamma_x);

    //--- Whiten the data using gamma_x and at the same time compute the
    //    logarithm of the determinant. On entry, A and y are equal to H and x.
    mactrick(gamma_x,lndeterminant,A,y);

    //-- Save y in vector r
    cblas_dcopy(m,y,1,r,1);

    //--- Simulate dgels using ATLAS subroutines (ordinary least-squares)
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,n,m,1.0,A,m,0.0,C_theta,n);
    clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);
    clapack_dpotri(CblasColMajor,CblasUpper,n,C_theta,n);
    cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A,m,y,1,0.0,dummyx,1);
    cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummyx,1,0.0,theta,1);

    //--- Compute y_hat = A*theta and add that to -1*y(=r) to get residuals
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        m,1,n,1.0,A,m,theta,n,-1.0,r,m);

    product = cblas_ddot(m,r,1,r,1);
    //--- compute sigma_eta
    sigma_eta = sqrt(product/static_cast<double>(m)); 

    //--- sigma_eta and lndeterminant are now set, copy their values
    //    to the local variables sigma_eta_ and lndeterminant_ which
    //    are returned to the caller.
    sigma_eta_     = sigma_eta;
    lndeterminant_ = lndeterminant;
  }
