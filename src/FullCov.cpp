/*! \file   FullCov.cpp
 *  \author Machiel Bos
 *
 * A class that computes the likelihood function, including the Least-Squares
 * part, using the standard method of using the full covariance matrix and
 * afterwards eliminating all the rows and columns for the missing data.
 *
 * \date 7/10/2012  Tomar
 */
//==============================================================================
  #include "FullCov.h"
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstdio>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!-----------------
  FullCov::FullCov(void)
//---!!-----------------
  {
    int           i,j;
    NoiseModel    *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- get number of noise model parameters to estimate
    Nparam = noisemodel->get_Nparam();
    N      = m-Ngaps;
#ifdef DEBUG
    cout << "Nparam=" << Nparam << ", m=" << m << ", n=" << n 
         << ", Ngaps=" << Ngaps << ", N=" << N << endl;
#endif

    //--- From the base class Likelihood, the values of m, n & Ngaps are known.
    try {
      gamma_x = new double[m];
      r       = new double[N];
      y       = new double[N];
      A       = new double[N*n];
      dummyH  = new double[N*n];
      dummyx  = new double[N];
      dummy   = new double[N];
      C       = new double[N*N];

      //--- Just avoid trouble
      memset(C,0.0,N*N*sizeof(double));
    }
    catch (bad_alloc) {
      cerr << "FullCov: Need more memory for variables!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Already eliminate rows of missing data from x and H and store the
    //    results in dummy_x and dummy_H.
    j=0;
    for (i=0;i<m;i++) {
      if (!isnan(x[i])) {
        dummyx[j] = x[i];
        cblas_dcopy(n,&H[i],m,&dummyH[j],N);
        j++;
      }
    }
    if (j!=N) {
      cerr << "Huh, something's wrong here: j=" << j << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Free up memory
 */
//---!!------------------
  FullCov::~FullCov(void)
//---!!------------------
  {
    if (gamma_x!=NULL)   delete[]  gamma_x;
    if (r!=NULL)         delete[]  r;
    if (y!=NULL)         delete[]  y;
    if (A!=NULL)         delete[]  A;
    if (dummyH!=NULL)    delete[]  dummyH;
    if (dummyx!=NULL)    delete[]  dummyx;
    if (dummy!=NULL)     delete[]  dummy;
    if (C!=NULL)         delete[]  C;
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
  void FullCov::compute_LeastSquares(double *param, double& 
					lndeterminant_, double& sigma_eta_)
//-------------------------------------------------------------------------
  {
    int         i,j,i0,j0,i1,j1,*ipiv;
    double      product,lambda,lambda_min,lambda_max;
    NoiseModel  *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- I assume the noise parameters, and thus the covariance matrix, have 
    //    been changed. As a result, matrix A and y need again to be computed.
    cblas_dcopy(N*n,dummyH,1,A,1);
    cblas_dcopy(N,dummyx,1,y,1);

    //--- Create pivots
    ipiv = new int[N];
    for (i=0;i<N;i++) ipiv[i]=i;

    //--- Create the Covariance matrix. 
    noisemodel->get_covariance(param,m,gamma_x);
#ifdef DEBUG
    show_matrix("gamma_x",gamma_x,m,1);
#endif

    i1=0;
    for (i0=0;i0<m;i0++) {
      if (!isnan(x[i0])) {
        j1=0;
        for (j0=0;j0<m-i0;j0++) {
          if (!isnan(x[i0+j0])) {
            if (i1>=N || j1>=N || (i1+(i1+j1)*N)>=N*N) {
              cerr << "Oops, i1=" << i1 << ", j1=" << j1 << endl;
              exit(EXIT_FAILURE);
            }
            C[i1+(i1+j1)*N] = gamma_x[j0];
            j1++;
          }
        }
        i1++;
      }
    }
    if (i1!=N) {
     cerr << "Something's wrong here: i1=" << i1 << endl;
     exit(EXIT_FAILURE);
    }
#ifdef DEBUG
    show_matrix("C",C,N,N);
#endif

    //--- Perform Cholesky decomposition, compute determinant and compute
    //    matrix A and vector y.
    clapack_dpotrf(CblasColMajor,CblasUpper,N,C,N);
    clapack_dgetrs(CblasColMajor,CblasTrans,N,n,C,N,ipiv,A,N);
    clapack_dgetrs(CblasColMajor,CblasTrans,N,1,C,N,ipiv,y,N);
#ifdef DEBUG
    show_matrix("A",A,N,n);
    show_matrix("y",y,N,1);
#endif

    //--- Compute logaritm of determinant
    lndeterminant = 0.0;
    for (i=0;i<N;i++) {
      lndeterminant += 2.0*log(C[i+N*i]);
    }

    //-- Save y in vector r
    cblas_dcopy(N,y,1,r,1);

    //--- Simulate dgels using ATLAS subroutines (ordinary least-squares)
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,n,N,1.0,A,N,0.0,C_theta,n);
    clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);
    clapack_dpotri(CblasColMajor,CblasUpper,n,C_theta,n);
    cblas_dgemv(CblasColMajor,CblasTrans,N,n,1.0,A,N,y,1,0.0,dummy,1);
    cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummy,1,0.0,theta,1);
#ifdef DEBUG
    show_matrix("theta",theta,n,1);
    show_matrix("C_theta",C_theta,n,n);
#endif

    //--- Compute y_hat = A*theta and add that to -1*y(=r) to get residuals
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        N,1,n,1.0,A,N,theta,n,-1.0,r,N);

    product = cblas_ddot(N,r,1,r,1);
    //--- compute sigma_eta
    sigma_eta = sqrt(product/static_cast<double>(N)); 

    //--- sigma_eta and lndeterminant are now set, copy their values
    //    to the local variables sigma_eta_ and lndeterminant_ which
    //    are returned to the caller.
    sigma_eta_     = sigma_eta;
    lndeterminant_ = lndeterminant;
#ifdef DEBUG
    cout << "sigma_eta     = " << sigma_eta << endl;
    cout << "product       = " << product << endl;
    cout << "lndeterminant = " << lndeterminant << endl;
#endif

    //--- free memory
    delete[] ipiv;
  }
