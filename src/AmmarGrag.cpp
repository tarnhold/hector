/*! \file   AmmarGrag.cpp
 *  \author Machiel Bos
 *
 * My implementation of the Fast Toeplitz solver of Ammar and Gragg (1988).
 *
 * References:
 * Ammar GS, Gragg WB (1988) Superfast solution of real positive definite 
 * Toeplitz systems. SIAM J Matrix Anal Appl, 9:61–76.

 * \date 9/2/2012  CIIMAR, Porto
 */
//==============================================================================
  #include "AmmarGrag.h"
  #include <iostream>
  #include <ostream>
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <sys/time.h>

//  #define DEBUG
//  #define TIME

//==============================================================================
// Subroutines
//==============================================================================


//---!!---------------------
  AmmarGrag::AmmarGrag(void)
//---!!---------------------
  {
    int           i,k,Status,nyc,offset;

    using namespace std;
    //--- Make use of the fact that FFT can be performed in parallel
    #if OMP == 1
      Nthreads = omp_get_max_threads();
    #else
      Nthreads = 1;
    #endif
    cout << endl << "Number of CPU's used (threads) = " << Nthreads 
         << endl << endl;

    //--- To speed things up, create permanent variable space
    ny = 2*m-1;
    nyc= ny/2+1;

    //--- For the more advanced ways of treating gaps in the data I need
    //    to have zero entries on the rows when there is a gap.
    if (Ngaps>0) {
      for (i=0;i<m;i++) {
        if (isnan(x[i])) {
          x[i] = 0.0;
          cblas_dscal(n,0.0,&H[i],m);
        }
      }
    }

    try {
      //--- For the least-squares we need the following variables
      gamma_x = new double[m];
      A1      = new double[m*n];
      A2      = new double[m*n];
      y1      = new double[m];
      y2      = new double[m];
      F_H     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc*n);
      F_x     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc*1);
      if (Ngaps>0) {
        F_F =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc*Ngaps);
        G1  = new double[m*Ngaps];
        G2  = new double[m*Ngaps];
        QA  = new double[Ngaps*n];
        Qy  = new double[Ngaps*1];
        Qt  = new double[Ngaps*1];
        M   = new double[Ngaps*Ngaps];
      } else {
        G1 = G2 = NULL;
        QA = Qy = Qt = M = NULL;
      }
      //--- For step 2 we need m+(m-1) variables. The FFT transform
      //    of real data requires an array size of m/2+1 
      l1      = new double[ny];
      l2      = new double[ny];
      F_l1    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc);
      F_l2    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc);
      dummy   = new double[ny*Nthreads];
      F_dummy = (fftw_complex*) 
			fftw_malloc(sizeof(fftw_complex) * nyc*Nthreads);
    }
    catch (bad_alloc()) {
      cerr << "AmmarGrag: need more memory!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- For step 2 we need to create a forward and backward plan
    plan_forward  = fftw_plan_dft_r2c_1d ( ny,   x, F_x, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d ( ny, F_x,   x, FFTW_ESTIMATE);

    //--- Already compute FFT of H, x and F -> F_H, F_x and F_F
    //    The matrices H, x and F are known through the base class Likelihood
    for (i=0;i<Nthreads;i++) {
      memset(&dummy[i*ny+m],0,(m-1)*sizeof(double)); //--- zero padding 
    }
    #pragma omp parallel for private(i,offset)
    for (i=0;i<n;i++) {
      #if OMP == 1
        offset = ny*(omp_get_thread_num ());
      #else 
        offset = 0; 
      #endif
      cblas_dcopy(m,&H[m*i],1,&dummy[offset],1);
      fftw_execute_dft_r2c(plan_forward,&dummy[offset],&F_H[nyc*i]);
    }
    cblas_dcopy(m,x,1,dummy,1);
    fftw_execute_dft_r2c(plan_forward,dummy,F_x);
    #pragma omp parallel for private(i,offset)
    for (i=0;i<Ngaps;i++) {
      #if OMP == 1
        offset = ny*(omp_get_thread_num ());
      #else
        offset = 0;
      #endif
      cblas_dcopy(m,&F[m*i],1,&dummy[offset],1);
      fftw_execute_dft_r2c(plan_forward,&dummy[offset],&F_F[nyc*i]);
    }
  }



/*! Free memory
 */
//---!!----------------------
  AmmarGrag::~AmmarGrag(void)
//---!!----------------------
  {
    int     i;

    if (gamma_x!=NULL)  delete[] A1;
    if (A1!=NULL)       delete[] A1;
    if (A2!=NULL)       delete[] A2;
    if (y1!=NULL)       delete[] y1;
    if (y2!=NULL)       delete[] y2;
    if (G1!=NULL)       delete[] G1;
    if (G2!=NULL)       delete[] G2;
    if (l1!=NULL)       delete[] l1;
    if (l2!=NULL)       delete[] l2;
    if (dummy!=NULL)    delete[] dummy;
    if (QA!=NULL)       delete[] QA;
    if (Qy!=NULL)       delete[] Qy;
    if (Qt!=NULL)       delete[] Qt; 

    //--- FFTW stuff
    fftw_destroy_plan ( plan_forward );
    fftw_destroy_plan ( plan_backward );
    fftw_free (F_H);
    fftw_free (F_x);
    fftw_free (F_F);
    fftw_free (F_l1);
    fftw_free (F_l2);
    fftw_free (F_dummy);
  }



/*! For the moment I use the slow AmmarGrag algorithm. To avoid too much
 *  copying around, I use l1 & l2 also as workspace. A side effect
 *  is that in the end, the pointers of l1 & l2 are switched when m is odd.
 *
 * This algorithm is based on chapter 3 of Ng's book on Toeplitz solvers
 * To be precise, page 28, the Durbin-Levinson algorithm
 *
 *  \param[in]  gamma_x   first column of covariance matrix
 *  \param[out] l1        first part of the Gohberg Semencul matrix
 *  \param[out] l2        second part of the Gohberg Semencul matrix
 */
//-------------------------------------------------------------------
  double AmmarGrag::step1(double *gamma_x, double **l1, double **l2)
//-------------------------------------------------------------------
  {
    int      i,j,Status;
    double   delta,ln_determinant,*dummy=NULL,gamma;
 
    using namespace std;
    //--- use r1 and r2 as working arrays (Note that [m..2m-1] is the
    //    needed zero padding. [0..2*m-1] is just to avoid random NaN's 
    memset(*l1,0,ny*sizeof(double));
    memset(*l2,0,ny*sizeof(double));

    //--- Start the algorithm
    delta = gamma_x[0];
    ln_determinant = log(delta);
    for (i=0;i<m-1;i++) {
      if (i==0) {
        gamma = -gamma_x[i+1]/delta;
      } else {
        gamma = -(gamma_x[i+1] + cblas_ddot(i,&gamma_x[1],1,
							&(*l2)[1],1))/delta;
        cblas_dcopy(i,&(*l2)[1],1,&(*l1)[2],1);
        cblas_daxpy(i,gamma,&(*l2)[1],-1,&(*l1)[2],1);

        //--- prepare next round
        dummy = *l2;
        *l2 = *l1;
        *l1 = dummy;
      }
      (*l2)[1] = gamma;
      delta = gamma_x[0] + cblas_ddot(i+1,&gamma_x[1],1,&(*l2)[1],-1);
      ln_determinant += log(delta);
    }

    //--- Create l1
    cblas_dcopy(m-1,&(*l2)[1],1,&(*l1)[1],-1);
    (*l1)[0] = 1.0;

    //--- Incorporate 1/delta_m into matrices l1 and l2.
    //    l1[m:2m-1] and l2[m:2m-1] are zero and don't need to be scaled.
    cblas_dscal(m,1.0/(static_cast<double>(1.0)*sqrt(delta)),*l1,1);
    cblas_dscal(m,1.0/(static_cast<double>(1.0)*sqrt(delta)),*l2,1);

    //--- Finally, already perform forward FFT's
    fftw_execute_dft_r2c(plan_forward,*l1,F_l1);
    fftw_execute_dft_r2c(plan_forward,*l2,F_l2);
    
    return ln_determinant;
  }



/*! multiply matrixl1 & l2 with matrix B. In the frequency domain this
 *  becomes equal to elementwise multiplication and inv(FFT).
 *
 *  B stands for A, y and G (design, observation and F-matrix)
 *
 * \param[in]   n_columns   number of columns of FFT-ed matrix F_B
 * \param[in]   F_B         FFT-ed matrix B
 * \param[out]  B1          result of multiplication of B*L1 in FFT-space
 * \param[out]  B2          result of multiplication of B*L2 in FFT-space
 */
//-------------------------------------------------------------------------
  void AmmarGrag::step2(int n_columns, fftw_complex *F_B, 
						    double *B1, double *B2)
//-------------------------------------------------------------------------
  { 
    int       i,j,offset1,offset2,nyc;
    double    Scale;

    using namespace std;
    //--- compute number of rows in the Fourier transformed variables
    Scale  = 1.0/static_cast<double>(ny);
    nyc = ny/2 + 1;

    //--- Perform FFT-multiplication for each column of A seperately
    #pragma omp parallel for private(i,j,offset1,offset2)
    for (i=0;i<n_columns;i++) {

      //--- Which thread is computing this?
      #if OMP == 1
        offset1 = nyc*(omp_get_thread_num ());
        offset2 = ny*(omp_get_thread_num ());
      #else
        offset1 = offset2 = 0;
      #endif

      for (j=0;j<nyc;j++) {
        F_dummy[offset1+j][0] = 
			F_l1[j][0]*F_B[nyc*i+j][0] - F_l1[j][1]*F_B[nyc*i+j][1];
        F_dummy[offset1+j][1] = 
			F_l1[j][0]*F_B[nyc*i+j][1] + F_l1[j][1]*F_B[nyc*i+j][0];
      }
      fftw_execute_dft_c2r(plan_backward,&F_dummy[offset1],&dummy[offset2]);
      cblas_dcopy(m,&dummy[offset2],1,&B1[m*i],1);
      cblas_dscal(m,Scale,&B1[m*i],1);

      for (j=0;j<nyc;j++) {
        F_dummy[offset1+j][0] = 
			F_l2[j][0]*F_B[nyc*i+j][0] - F_l2[j][1]*F_B[nyc*i+j][1];
        F_dummy[offset1+j][1] = 
			F_l2[j][0]*F_B[nyc*i+j][1] + F_l2[j][1]*F_B[nyc*i+j][0];
      }
      fftw_execute_dft_c2r(plan_backward,&F_dummy[offset1],&dummy[offset2]);
      cblas_dcopy(m,&dummy[offset2],1,&B2[m*i],1);
      cblas_dscal(m,Scale,&B2[m*i],1);
    }
  }



/*! compute least-squares using given values of noise parameters.
 *
 * \param[in]  param           arrray with noise parameters
 * \param[out] lndeterminant_  logarithm of likelihood value
 * \param[out] sigma_eta_      standard deviation of the driving white noise
 */
//-----------------------------------------------------------------------
  void AmmarGrag::compute_LeastSquares(double *param,
                                double& lndeterminant, double& sigma_eta)
//-----------------------------------------------------------------------
  {
    int            i,j,*ipiv;
    double         product;
    NoiseModel     *noisemodel=NoiseModel::getInstance();
    struct timeval start,end,start_total;
    long           mtime, seconds, useconds;  


    using namespace std;
    gettimeofday(&start_total, NULL);

    //--- Create pivots
    ipiv = new int[Ngaps];
    for (i=0;i<Ngaps;i++) ipiv[i]=i;

    //--- Create the 1st column of the Covariance matrix. Since this is
    //    a Toeplitz matrix, this column vector is sufficient.
    //    Note that this is for unit variance value of sigma_eta.
    noisemodel->get_covariance(param,m,gamma_x);

    //--- Perform step 1 : compute l1 and l2 vectors (including delta_m)
    gettimeofday(&start, NULL);
    lndeterminant = step1(gamma_x,&l1,&l2);
    gettimeofday(&end, NULL);
#ifdef TIME
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "step1: " << mtime << " msec" << endl;
#endif

    //--- Compute the auxiliary matrices A1, A2, y1, y2, G1 and G2;
    gettimeofday(&start, NULL);
    step2(n,F_H,A1,A2);
    step2(1,F_x,y1,y2);
    if (Ngaps>0) step2(Ngaps,F_F,G1,G2);
    gettimeofday(&end, NULL);
#ifdef TIME
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "step2: " << mtime << " msec" << endl;
#endif

    if (Ngaps==0) {
      //--- Compute C_theta
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						n,m,1.0,A2,m, 0.0,C_theta,n);
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						n,m,1.0,A1,m,-1.0,C_theta,n);

      //--- Now compute inv(A1'*A1 - A2'*A2)
      clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);
      clapack_dpotri(CblasColMajor,CblasUpper,n,C_theta,n);

      //--- Compute A1'*y1 - A2'*y2
      cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A2,m,y2,1, 0.0,dummy,1);
      cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A1,m,y1,1,-1.0,dummy,1);
      //--- Finally, perform Least-Squares
      cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummy,1,0.0,theta,1);

      //--- Compute residuals
      cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A1,m,theta,1,-1.0,y1,1);
      cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A2,m,theta,1,-1.0,y2,1);
    
      product = cblas_ddot(m,y1,1,y1,1) - cblas_ddot(m,y2,1,y2,1);

      //------------ End No-Gap case ------------------------

    } else {
      //--- Compute M = Chol(G'*G)
      memset(M,0,(Ngaps*Ngaps)*sizeof(double)); // put M to zero
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						Ngaps,m,1.0,G2,m, 0.0,M,Ngaps);
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						Ngaps,m,1.0,G1,m,-1.0,M,Ngaps);
      clapack_dpotrf(CblasColMajor,CblasUpper,Ngaps,M,Ngaps);

      //--- Adjust lndeterminant
#ifdef DEBUG
      cout << "lndeterminant=" << lndeterminant << endl;
#endif
      for (i=0;i<Ngaps;i++) {
        lndeterminant += 2.0*log(M[i+i*Ngaps]);
      }
#ifdef DEBUG
      cout << " after correction: lndeterminant=" << lndeterminant << endl;
#endif

      //--- Solve M'*Q_A=G'*A
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,n,m,1.0,G2,m,A2,m, 0.0,QA,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,n,m,1.0,G1,m,A1,m,-1.0,QA,Ngaps);
      clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,n,M,Ngaps,ipiv,QA,Ngaps);
      //--- Solve M'*Q_y=G'*y
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G2,m,y2,m, 0.0,Qy,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G1,m,y1,m,-1.0,Qy,Ngaps);
      clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,1,M,Ngaps,ipiv,Qy,Ngaps);

      //--- Compute C_theta
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,Ngaps,1.0,QA,Ngaps,0.0,C_theta,n);
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,m, 1.0,A2,m,1.0,C_theta,n);
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,m,1.0,A1,m,-1.0,C_theta,n);

      //--- Now compute inv(A'*A)
      clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);
      clapack_dpotri(CblasColMajor,CblasUpper,n,C_theta,n);

      //--- Compute A'*y - Q_A'*Q_y
      cblas_dgemv(CblasColMajor,CblasTrans,
				Ngaps,n,1.0,QA,Ngaps,Qy,1, 0.0,dummy,1);
      cblas_dgemv(CblasColMajor,CblasTrans,
				m,n,1.0,A2,m,y2,1, 1.0,dummy,1);
      cblas_dgemv(CblasColMajor,CblasTrans,
				m,n,1.0,A1,m,y1,1,-1.0,dummy,1);

      //--- Finally, perform Least-Squares
      cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummy,1,0.0,theta,1);

      //--- Compute residuals t1=A1*theta-y1 -> stored in y1
      cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A1,m,theta,1,-1.0,y1,1);
      cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A2,m,theta,1,-1.0,y2,1);
      //--- Solve M'*Q_t=G'*(A*theta-y)
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G2,m,y2,m, 0.0,Qt,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G1,m,y1,m,-1.0,Qt,Ngaps);
      clapack_dgetrs(CblasColMajor,CblasTrans,Ngaps,1,M,Ngaps,ipiv,Qt,Ngaps);

      product = cblas_ddot(m,y1,1,y1,1) - cblas_ddot(m,y2,1,y2,1) -
						cblas_ddot(Ngaps,Qt,1,Qt,1);

    }

    //--- compute sigma_eta
    sigma_eta = sqrt(product/static_cast<double>(m-Ngaps));

    //--- free memory
    delete[] ipiv;

#ifdef DEBUG
    cout << "product=" << product << endl;
    cout << "sigma_eta =" << sigma_eta << endl;
#endif
#ifdef TIME
    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start_total.tv_sec;
    useconds = end.tv_usec - start_total.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "Total LS: " << mtime << " msec" << endl;
#endif

  }
