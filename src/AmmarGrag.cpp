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
    int           i,j,nyc,offset;

    using namespace std;
    //--- Make use of the fact that FFT can be performed in parallel
    Nthreads = omp_get_max_threads();
    if (Nthreads>8) {
      omp_set_dynamic(0);     // Explicitly disable dynamic teams
      omp_set_num_threads(8); // Use no more than 8 threads for all OMP stuff
      Nthreads = omp_get_max_threads();
    }
    cout << endl << "Number of CPU's used (threads) = " << Nthreads 
         << endl << endl;

    //--- To speed things up, create permanent variable space
    ny = 2*m-1;
    nyc= ny/2+1;

    //--- For the more advanced ways of treating gaps in the data I need
    //    to have zero entries on the rows when there is a gap.
    if (Ngaps>0) {
      for (i=0;i<m;i++) {
        if (std::isnan(x[i])) {
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
        G1  = new double[m*Ngaps];
        G2  = new double[m*Ngaps];
        QA  = new double[Ngaps*n];
        Qy  = new double[Ngaps*1];
        Qt  = new double[Ngaps*1];
        M   = new double[Ngaps*Ngaps];
        M2  = new double[Ngaps*Ngaps];
      } else {
        G1 = G2 = NULL;
        QA = Qy = Qt = M = M2 = NULL;
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

    //--- Already compute FFT of H, x -> F_H and F_x
    //    The matrices H, x are known through the base class MLEBase
    for (i=0;i<Nthreads;i++) {
      memset(&dummy[i*ny+m],0,(m-1)*sizeof(double)); //--- zero padding 
    }
    #pragma omp parallel for private(i,offset)
    for (i=0;i<n;i++) {
      offset = ny*(omp_get_thread_num ());
      cblas_dcopy(m,&H[m*i],1,&dummy[offset],1);
      fftw_execute_dft_r2c(plan_forward,&dummy[offset],&F_H[nyc*i]);
    }
    cblas_dcopy(m,x,1,dummy,1);
    fftw_execute_dft_r2c(plan_forward,dummy,F_x);
    //--- Trying to make use of the many zeros in F. Array index contains
    //    for each column the row index where the first 1 is encountered.
    if (Ngaps>0) {
      index = new int[Ngaps];
      for (i=0;i<Ngaps;i++) {
        j=0;
        while (F[j + m*i]<0.99) j++;
        if (j>=m) {
          cerr << "Trouble in paradise! j=" << j << ", m=" << m << endl;
          exit(EXIT_FAILURE);
        }
        index[i] = j;
      }
    } else {
      index = NULL;
    }
  }



/*! Free memory
 */
//---!!----------------------
  AmmarGrag::~AmmarGrag(void)
//---!!----------------------
  {
    if (gamma_x!=NULL)  delete[] gamma_x;
    if (A1!=NULL)       delete[] A1;
    if (A2!=NULL)       delete[] A2;
    if (y1!=NULL)       delete[] y1;
    if (y2!=NULL)       delete[] y2;
    if (G1!=NULL)       delete[] G1;
    if (G2!=NULL)       delete[] G2;
    if (M!=NULL)        delete[] M;
    if (M2!=NULL)       delete[] M2;
    if (l1!=NULL)       delete[] l1;
    if (l2!=NULL)       delete[] l2;
    if (dummy!=NULL)    delete[] dummy;
    if (QA!=NULL)       delete[] QA;
    if (Qy!=NULL)       delete[] Qy;
    if (Qt!=NULL)       delete[] Qt;
    if (index!=NULL)    delete[] index; 

    //--- FFTW stuff
    fftw_destroy_plan ( plan_forward );
    fftw_destroy_plan ( plan_backward );
    fftw_free (F_H);
    fftw_free (F_x);
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
 *  \return logarithm of determinant of Covariance matrix
 */
//-------------------------------------------------------------------
  double AmmarGrag::step1(double *gamma_x, double **l1, double **l2)
//-------------------------------------------------------------------
  {
    int      i;
    double   delta,ln_determinant_C,*dummy=NULL,gamma;
 
    using namespace std;
    //--- use r1 and r2 as working arrays (Note that [m..2m-1] is the
    //    needed zero padding. [0..2*m-1] is just to avoid random NaN's 
    memset(*l1,0,ny*sizeof(double));
    memset(*l2,0,ny*sizeof(double));

    //--- Start the algorithm
    delta = gamma_x[0];
    ln_determinant_C = log(delta);
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
      ln_determinant_C += log(delta);
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
    
    return ln_determinant_C;
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
      offset1 = nyc*(omp_get_thread_num ());
      offset2 = ny*(omp_get_thread_num ());

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



/*! compute l1 and l2. At the same time already compute G1, G2, y1 and y2.
 *
 * \param[in]  param             arrray with noise parameters
 * \param[out] ln_determinant_C  logarithm of determinant of Covariance matrix
 */
//-------------------------------------------------
  void AmmarGrag::prepare_covariance(double *param)
//-------------------------------------------------
  {
    NoiseModel     &noisemodel=NoiseModel::getInstance();
    struct timeval start,end;
    int            i;
#ifdef TIME
    long           mtime, seconds, useconds;
#endif

    using namespace std;
    //--- Create the 1st column of the Covariance matrix. Since this is
    //    a Toeplitz matrix, this column vector is sufficient.
    //    Note that this is for unit variance value of sigma_eta.
    noisemodel.get_covariance(param,m,gamma_x);

    //--- Perform step 1 : compute l1 and l2 vectors (including delta_m)
    gettimeofday(&start, NULL);
    ln_det_C = step1(gamma_x,&l1,&l2); //-- A MLEBase-class variable
    gettimeofday(&end, NULL);
#ifdef TIME
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "step1: " << mtime << " msec" << endl;
#endif

    gettimeofday(&start, NULL);
    //--- It's faster to copy directly the result into G1 and G2 then to
    //    perform FFT and iFFT's.
    if (Ngaps>0) {
      memset(G1,0,(m*Ngaps)*sizeof(double));
      memset(G2,0,(m*Ngaps)*sizeof(double));
      #pragma omp parallel for private(i)
      for (i=0;i<Ngaps;i++) {
        cblas_dcopy(m-index[i],l1,1,&G1[index[i]+i*m],1);
        cblas_dcopy(m-index[i],l2,1,&G2[index[i]+i*m],1);
      }
    }
    
    gettimeofday(&end, NULL);
#ifdef TIME
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "step2: " << mtime << " msec" << endl;
#endif

    if (Ngaps>0) {
      //--- Compute M = Chol(G'*G)
      memset(M ,0,(Ngaps*Ngaps)*sizeof(double)); // put M to zero
      memset(M2,0,(Ngaps*Ngaps)*sizeof(double)); // put M2 to zero
      #pragma omp parallel sections
      {
        #pragma omp section
        cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						Ngaps,m,1.0,G2,m,0.0,M2,Ngaps);
        #pragma omp section
        cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						Ngaps,m,1.0,G1,m,0.0,M,Ngaps);
      }
      cblas_daxpy(Ngaps*Ngaps,-1.0,M2,1,M,1);

      //--- Perform Cholesky decomposition of M
      clapack_dpotrf(CblasColMajor,CblasUpper,Ngaps,M,Ngaps);

      //--- Adjust ln_det_C which is a MLEBase-class variable
#ifdef DEBUG
      cout << "ln_det_C=" << ln_det_C << endl;
#endif
      for (i=0;i<Ngaps;i++) {
        ln_det_C += 2.0*log(M[i+i*Ngaps]);
      }
#ifdef DEBUG
      cout << " after correction: ln_det_C="<<ln_det_C << endl;
#endif

    }
  }



/*! compute least-squares using given values of noise parameters.
 *
 * \param[in]  param             arrray with noise parameters
 * \param[out] ln_determinant_I  logarithm of determinant of information matrix
 * \param[out] sigma_eta         standard deviation of the driving white noise
 */
//---------------------------------------------------
  void AmmarGrag::compute_LeastSquares(double * /*param*/)
//---------------------------------------------------
  {
    int            i,*ipiv,k;
    double         product;

    using namespace std;
    //--- Create pivots (check which is larger n or Ngaps)
    k = n;
    if (Ngaps>n) k=Ngaps;
    ipiv = new int[k];
    for (i=0;i<k;i++) ipiv[i]=i;

    //--- Compute the auxiliary matrices A1, A2. These are always recomputed
    //    because I change H when looking for offsets.
    step2(1,F_x,y1,y2);
    step2(n,F_H,A1,A2);

    if (Ngaps==0) {
      //--- Compute C_theta
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						n,m,1.0,A2,m, 0.0,C_theta,n);
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
						n,m,1.0,A1,m,-1.0,C_theta,n);

      //--- Now compute inv(A1'*A1 - A2'*A2)
      clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);

      //--- Now we can compute inverse of C_theta.  
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

      //--- Now compute inv(A'*A). First compute Cholesky decomposition
      clapack_dpotrf(CblasColMajor,CblasUpper,n,C_theta,n);

      //--- Now we can compute inverse of C_theta.  
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

    //--- compute sigma_eta which is a MLEBase-class variable
    sigma_eta = sqrt(product/static_cast<double>(m-Ngaps));

    //--- free memory
    delete[] ipiv;

#ifdef DEBUG
    cout << "product=" << product << endl;
    cout << "sigma_eta =" << sigma_eta << endl;
#endif
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
