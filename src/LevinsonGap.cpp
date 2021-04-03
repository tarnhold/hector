/*! \file   LevinsonGap.cpp
 *  \author Machiel Bos
 *
 * A class that computes the likelihood function, including the Least-Squares
 * part, using the Levinson method to perform a Cholesky decomposition of the
 * Covariance matric. This is my new approach, assuming that there are
 * gaps.
 *
 * \date 3/2/2012  CIIMAR, Porto
 */
//==============================================================================
  #include "LevinsonGap.h"
  #include <cmath>
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstdio>
  #include <ctime>

//  #define DEBUG
//  #define TIME

//==============================================================================
// Subroutines
//==============================================================================



//---!!-------------------------
  LevinsonGap::LevinsonGap(void)
//---!!-------------------------
  {
    char          NoTrans='N';
    int           i;
    NoiseModel    *noisemodel=NoiseModel::getInstance();

    using namespace std;
    //--- get number of noise model parameters to estimate
    Nparam = noisemodel->get_Nparam();

    //--- From the base class Likelihood, the values of n & m are known.
    try {
      gamma_x = new double[m];
      r       = new double[m];
      y       = new double[m];
      A       = new double[m*n];
      G       = new double[m*Ngaps];
      t       = new double[m];
      dummyH  = new double[m*n];
      dummyx  = new double[m];
      dummyF  = new double[m*Ngaps];
      U       = new double[2*m];
      V       = new double[2*m];
      Qy      = new double[Ngaps*1];
      QA      = new double[Ngaps*n];
      Qt      = new double[Ngaps*1];
      M       = new double[Ngaps*Ngaps];
      if (n>Ngaps) {
        dummy   = new double[m*n];
      } else {
        dummy   = new double[m*Ngaps];
      }
    }
    catch (bad_alloc) {
      cerr << "Levinson: Need more memory for variables!" << endl;
      exit(EXIT_FAILURE);
    }
    //--- Avoid trouble with random NaN's in memory
    memset(U,0.0,(2*m)*sizeof(double)); 
    memset(V,0.0,(2*m)*sizeof(double)); 
    memset(Qy,0.0,(Ngaps*1)*sizeof(double));
    memset(QA,0.0,(Ngaps*n)*sizeof(double));
    memset(Qt,0.0,(Ngaps*1)*sizeof(double));
    memset(Qt,0.0,(Ngaps*1)*sizeof(double));
    memset(M,0.0,(Ngaps*Ngaps)*sizeof(double)); 

    //--- For the more advanced ways of treating gaps in the data I need
    //    to have zero entries on the rows when there is a gap.
    if (Ngaps>0) {
      for (i=0;i<m;i++) {
        if (isnan(x[i])) {
          x[i] = 0.0;
          cblas_dscal(n,0.0,&H[i],m);
        }
      }
    } else {
      cerr << "LevisonGap:: Ngaps==0! wrong class called...." << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Free up memory
 */
//---!!--------------------------
  LevinsonGap::~LevinsonGap(void)
//---!!--------------------------
  {
    if (gamma_x!=NULL)   delete[]  gamma_x;
    if (r!=NULL)         delete[]  r;
    if (y!=NULL)         delete[]  y;
    if (A!=NULL)         delete[]  A;
    if (G!=NULL)         delete[]  G;
    if (t!=NULL)         delete[]  t;
    if (M!=NULL)         delete[]  M;
    if (dummy!=NULL)     delete[]  dummy;
    if (dummyH!=NULL)    delete[]  dummyH;
    if (dummyx!=NULL)    delete[]  dummyx;
    if (dummyF!=NULL)    delete[]  dummyF;
    if (QA!=NULL)        delete[]  QA;
    if (Qy!=NULL)        delete[]  Qy;
    if (Qt!=NULL)        delete[]  Qt;
    if (U!=NULL)         delete[]  U;
    if (V!=NULL)         delete[]  V;
  }



/*! Mactrick
 *
 * \param[in]  gamma_x       : autocovariance matrix of noise
 * \param[out] lndeterminant : logarithm of determinant of C
 * \param[out] A             : U*H = A where U'*U=Cinv
 * \param[out] y             : U*x = y
 * \param[out] G             : U*F = G
 */
//---------------------------------------------------------------------------
  void LevinsonGap::mactrick(double *gamma_x,
                     double& lndeterminant_, double *A, double *y, double *G)
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
    cblas_dcopy(Ngaps,F,m,G,m);
    cblas_dscal(Ngaps,1.0/U[0],G,m);
    y[0] = x[0]/U[0];
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        m-1,n,1,1.0,U+1,m,A,m,0.0,dummyH+1,m);
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                        m-1,1,1,1.0,U+1,m,y,m,0.0,dummyx+1,m);
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                    m-1,Ngaps,1,1.0,U+1,m,G,m,0.0,dummyF+1,m);

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
      cblas_dcopy(Ngaps,F+k,m,G+k,m);
      cblas_daxpy(Ngaps,-1.0,dummyF+k,m,G+k,m);
      cblas_dscal(Ngaps,1.0/U[k + m*k_new],G+k,m);
      y[k] = (x[k]-dummyx[k])/U[k + m*k_new];
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                m-k-1,n,1,1.0,&U[k+1 + m*k_new],m,A+k,m,0.0,dummy+(k+1),m);
      cblas_daxpy(n*m,1.0,dummy,1,dummyH,1);
      cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                m-k-1,Ngaps,1,1.0,&U[k+1 + m*k_new],m,G+k,m,0.0,dummy+(k+1),m);
      cblas_daxpy(Ngaps*m,1.0,dummy,1,dummyF,1);
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
  void LevinsonGap::compute_LeastSquares(double *param, double& 
					lndeterminant_, double& sigma_eta_)
//-------------------------------------------------------------------------
  {
    int         i,j,*ipiv=NULL;
    double      product;
    NoiseModel  *noisemodel=NoiseModel::getInstance();
    time_t      start,end,start_total,end_total;


    using namespace std;
    time(&start_total);
    //--- Create pivots
    ipiv = new int[Ngaps];
    for (i=0;i<Ngaps;i++) ipiv[i]=i;

    //--- I assume the noise parameters, and thus the covariance matrix, have 
    //    been changed. As a result, matrix A and y need again to be computed.
    cblas_dcopy(m*n,H,1,A,1);
    cblas_dcopy(m,x,1,y,1);
    cblas_dcopy(m*Ngaps,F,1,G,1);

    //--- Create the 1st column of the Covariance matrix. Since this is
    //    a Toeplitz matrix, this column vector is sufficient.
    //    Note that this is for unit variance value of sigma_eta.
#ifdef DEBUG
    if (noisemodel->get_Nparam()>1) {
      cout << "param[0]=" << param[0] << ", param[1]=" << param[1] << endl;
    } else if (noisemodel->get_Nparam()>0) {
      cout << "param[0]=" << param[0] << endl;
    }
#endif
    noisemodel->get_covariance(param,m,gamma_x);

    //--- Whiten the data using gamma_x and at the same time compute the
    //    logarithm of the determinant. On entry, A and y are equal to H and x.
    //    G is equal to F.
    time(&start);
    mactrick(gamma_x,lndeterminant,A,y,G);
    time(&end);
#ifdef TIME
    cout << "mactrick: " << difftime(end,start) << " sec" << endl;
#endif
#ifdef DEBUG
    show_matrix("A",A,m,n);
    show_matrix("y",y,m,1);
    show_matrix("G",G,m,Ngaps);
#endif

    //--- Compute M = Chol(G'*G)
    time(&start);
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						Ngaps,m,1.0,G,m,0.0,M,Ngaps);
    clapack_dpotrf(CblasColMajor,CblasUpper,Ngaps,M,Ngaps);
    time(&end);
#ifdef TIME
    cout << "computation M: " << difftime(end,start) << " sec" << endl;
#endif
#ifdef DEBUG
    show_matrix("M",M,Ngaps,Ngaps);
#endif

    //--- Adjust lndeterminant
#ifdef DEBUG
    cout << "lndeterminant=" << lndeterminant << endl;
#endif
    for (i=0;i<Ngaps;i++) {
      lndeterminant += 2.0*log(M[i+i*Ngaps]);
    }
#ifdef DEBUG
    cout << "lndeterminant=" << lndeterminant << endl;
#endif

    //--- Solve M'*Q_A=G'*A
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,n,m,1.0,G,m,A,m,0.0,QA,Ngaps);
    clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,n,M,Ngaps,ipiv,QA,Ngaps);
    //--- Solve M'*Q_y=G'*y
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G,m,y,m,0.0,Qy,Ngaps);
    clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,1,M,Ngaps,ipiv,Qy,Ngaps);

#ifdef DEBUG
    show_matrix("QA",QA,Ngaps,n);
    show_matrix("Qy",Qy,Ngaps,1);
#endif

    //--- Compute C_theta
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,Ngaps,1.0,QA,Ngaps,0.0,C_theta,n);
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,m,1.0,A,m,-1.0,C_theta,n);
   
    //--- Now compute inv(A'*A)
    clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);
    clapack_dpotri(CblasColMajor,CblasUpper,n,C_theta,n);
#ifdef DEBUG
    show_matrix("C_theta",C_theta,n,n);
#endif
 
    //--- Compute A'*y - Q_A'*Q_y
    cblas_dgemv(CblasColMajor,CblasTrans,Ngaps,n,1.0,QA,Ngaps,Qy,1,0.0,dummy,1);
    cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A,m,y,1,-1.0,dummy,1);
#ifdef DEBUG
    show_matrix("dummy",dummy,n,1);
#endif

    //--- Finally, perform Least-Squares
    cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummy,1,0.0,theta,1); 
#ifdef DEBUG
    show_matrix("theta",theta,n,1);
#endif

    //--- Compute residuals
    cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A,m,theta,1,-1.0,y,1);
    //--- Solve M'*Q_t=G'*(A*theta-y)
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G,m,y,m,0.0,Qt,Ngaps);
    clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,1,M,Ngaps,ipiv,Qt,Ngaps);
    
#ifdef DEBUG
    show_matrix("t",y,m,1);
    show_matrix("Q_t",Qt,Ngaps,1);
#endif
    
    product = cblas_ddot(m,y,1,y,1) - cblas_ddot(Ngaps,Qt,1,Qt,1);
#ifdef DEBUG
    cout << "product=" << product << endl;
#endif

    //--- compute sigma_eta
    sigma_eta = sqrt(product/static_cast<double>(m-Ngaps)); 

    // Copy sigma_eta and lndeterminant to the variables sigma_eta_ and 
    // lndeterminant_ which are returned to the caller.
    sigma_eta_     = sigma_eta;
    lndeterminant_ = lndeterminant;

    //--- free memory
    if (Ngaps>0) delete[] ipiv;

    time(&end_total);
#ifdef TIME
    cout << "Total LS: " << difftime(end_total,start_total) << " sec" << endl;
#endif
  }
