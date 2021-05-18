/*! \file   AmmarGrag.cpp
 *  \author Machiel Bos
 *
 *  My implementation of the Fast Toeplitz solver of Ammar and Gragg (1988).
 *
 *  Implementation of the ARFIMA noise model using the tricks of Sowell (1992),
 *  Doornik and Ooms (2003) and Zinde-Wash (1988).
 *
 * References:
 * -----------
 *  Ammar GS, Gragg WB (1988) Superfast solution of real positive definite 
 *  Toeplitz systems. SIAM J Matrix Anal Appl, 9:61–76.
 *
 *  This script is part of Hector 1.9
 *
 *  Hector is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  Hector is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Hector.  If not, see <http://www.gnu.org/licenses/>
 *
 */
//==============================================================================
  #include "AmmarGrag.h"
  #include <iostream>
  #include <ostream>
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <algorithm>
  #include <sys/time.h>

//  #define DEBUG
//  #define TIME

//==============================================================================
// Subroutines
//==============================================================================


//---!!-------------------------------------------
  AmmarGrag::AmmarGrag(void) : tpi(8.0*atan(1.0)),
                               LARGE(9.9e99),
			       EPS(1.0e-6)
//---!!-------------------------------------------
  {
    int           i,j,k,Status,nyc,offset;

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
      t1      = new double[m];
      t2      = new double[m];
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
    if (t1!=NULL)       delete[] t1;
    if (t2!=NULL)       delete[] t2;
    if (G1!=NULL)       delete[] G1;
    if (G2!=NULL)       delete[] G2;
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
    int          i,j,Status;
    double       ln_determinant_C,*dummy=NULL,delta,gamma;
 
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
      if (std::isnan(delta) or fabs(delta)<1.0e-50) {
        cout << "Matrix does not appear to be positive definite..." << endl;

        cout << endl << "If you are using the GGM noise model, then read" <<
               " section 7.5 and 7.6 of the manual." << endl <<
               "You should increase the value of GGM_1mphi" << endl; 

        if (std::isnan(delta)==true) {
          cout << "Delta is NaN" << endl;
        } else {
          cout << "Delta is :" << delta << endl;
        }
        for (j=0;j<10;j++) {
          cout << "i=" << i << ", j=" << j << ", gamma_x=" << gamma_x[j] << 
		", l1=" << (*l1)[j] <<  ", l2=" << (*l2)[j] << endl;
        }
        exit(EXIT_FAILURE);
      }
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



/*! compute least-squares using given values of noise parameters.
 *
 * \param[in]  param             arrray with noise parameters
 * \param[out] ln_determinant_I  logarithm of determinant of information matrix
 * \param[out] sigma_eta         standard deviation of the driving white noise
 */
//---------------------------------------------------
  void AmmarGrag::compute_LeastSquares(double *param)
//---------------------------------------------------
  {
    NoiseModel   &noisemodel=NoiseModel::getInstance();
    int          i,j,*ipiv,k,info,one=1;
    char         Trans='T',Up='U';
    double       product,ms;
    clock_t      start,end,start2,end2;

    using namespace std;
    //--- Create pivots (check which is larger n or Ngaps)
    k = n;
    if (Ngaps>n) k=Ngaps;
    ipiv = new int[k];
    for (i=0;i<k;i++) ipiv[i]=i+1;

    //--- Create the 1st column of the Covariance matrix. Since this is
    //    a Toeplitz matrix, this column vector is sufficient.
    //    Note that this is for unit variance value of sigma_eta.
    noisemodel.get_covariance(param,m,gamma_x);

    //--- Perform step 1 : compute l1 and l2 vectors (including delta_m)
    ln_det_C = step1(gamma_x,&l1,&l2); //-- A MLEBase-class variable

    step2(1,F_x,y1,y2);
    step2(n,F_H,A1,A2);
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

    //--- Compute C_thetaInv = A1'*A1 - A2'*A2
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,m,1.0,A2,m, 0.0,C_thetaInv,n);
    cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,m,1.0,A1,m,-1.0,C_thetaInv,n);

    //--- Compute dummy = A1'*y1 - A2'*y2
    cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A2,m,y2,1, 0.0,dummy,1);
    cblas_dgemv(CblasColMajor,CblasTrans,m,n,1.0,A1,m,y1,1,-1.0,dummy,1);

    if (Ngaps>0) {
      //--- Compute matrix M
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
      dpotrf_(&Up,&Ngaps,M,&Ngaps,&info);

      //--- Adjust ln_det_C which is a MLEBase-class variable
      for (i=0;i<Ngaps;i++) {
        ln_det_C += 2.0*log(M[i+i*Ngaps]);
      }

      //--- Solve M'*Q_A=G1'*A1 - G2'*A2
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,n,m,1.0,G2,m,A2,m, 0.0,QA,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,n,m,1.0,G1,m,A1,m,-1.0,QA,Ngaps);
      dgetrs_(&Trans,&Ngaps,&n,M,&Ngaps,ipiv,QA,&Ngaps,&info);

      //--- Solve M'*Q_y=G1'*y1 - G2'*y2
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G2,m,y2,m, 0.0,Qy,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G1,m,y1,m,-1.0,Qy,Ngaps);
      dgetrs_(&Trans,&Ngaps,&one,M,&Ngaps,ipiv,Qy,&Ngaps,&info);

      //--- Adjust C_thetaInv with -QA'*QA
      cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,
					n,Ngaps,-1.0,QA,Ngaps,1.0,C_thetaInv,n);
      //--- Adjust dummy with -Q_A'*Q_y
      cblas_dgemv(CblasColMajor,CblasTrans,
				Ngaps,n,-1.0,QA,Ngaps,Qy,1,1.0,dummy,1);
    }

    //--- Now compute inv(A1'*A1 - A2'*A2 - [QA'*QA])
    cblas_dcopy(n*n,C_thetaInv,1,C_theta,1);
    dpotrf_(&Up,&n,C_theta,&n,&info);

    //--- Compute ln_det_I (The logarithm of the determinant of the
    //    Fisher Information matrix). det(C_theta^{-1}) = 1/det(C_theta)
    //    However, I = C_theta^{-1} so the inverses cancel.
    ln_det_I = 0.0;
    for (i=0;i<n;i++) {
      ln_det_I += 2.0*log(C_theta[i+i*n]);
    }
      
    //--- Now we can compute inverse of C_theta.  
    dpotri_(&Up,&n,C_theta,&n,&info);

    //--- Perform Least-Squares
    cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,C_theta,n,dummy,1,0.0,theta,1);

    //--- Compute residuals t1 and t2
    cblas_dcopy(m,y1,1,t1,1);
    cblas_dcopy(m,y2,1,t2,1);
    cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A1,m,theta,1,-1.0,t1,1);
    cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,A2,m,theta,1,-1.0,t2,1);
    
    product = cblas_ddot(m,t1,1,t1,1) - cblas_ddot(m,t2,1,t2,1);

    //--- For Gaps, subtract also Qt'*Qt 
    if (Ngaps>0) {
      //--- Solve M'*Q_t=(G1'*(A1*theta-y1) - G2'*(A2*theta-y2))
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G2,m,t2,m, 0.0,Qt,Ngaps);
      cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
					Ngaps,1,m,1.0,G1,m,t1,m,-1.0,Qt,Ngaps);
      dgetrs_(&Trans,&Ngaps,&one,M,&Ngaps,ipiv,Qt,&Ngaps,&info);
      product -= cblas_ddot(Ngaps,Qt,1,Qt,1);
    }

    //--- compute sigma_eta which is a MLEBase-class variable
    sigma_eta = sqrt(product/static_cast<double>(m-Ngaps));

    //--- Fisher Information matrix was computed using unit covariance
    //    matrix. Update with right scaling.
    ln_det_I -= 2.0*n*log(sigma_eta);

    //--- free memory
    delete[] ipiv;
  }



/*! Compute BIC_c for all possible offset locations. I assume that an offset
 *  was placed at the first observation (new column in H) and that
 *  compute_LeastSquares was called so that l1, l2 and a lot of other matrices
 *  and vectors already have been computed which are re-used here.
 */
//---------------------------------------------
  void AmmarGrag::compute_BIC_cs(double *BIC_c)
//---------------------------------------------
  {
    using namespace std;
    vector<double>             offsets;
    vector< vector<double> >   off_omp;
    NoiseModel      &noisemodel=NoiseModel::getInstance();
    int             i,j,l,column,ny,nyc,*ipiv,k,n_offsets,offset_index,N;
    int             one=1,info;
    char            Trans='T',Up='U';
    double          product,ms,time0,time1,lnf_s,lnf_theta;
    clock_t         start,end,start2,end2;
    bool            already_used;
    DesignMatrix    *designmatrix = DesignMatrix::getInstance();
    Observations    &observations = Observations::getInstance();
    double          *H_omp,*z_omp,*theta_omp,*C_theta_omp,*C_thetaInv_omp;
    double          *t1_omp,*t2_omp,*QA_omp,*Qt_omp,*dummy_omp,*A1_omp,*A2_omp;
    double          *t_ori,*x_ori,lambda,fact_k,sigma_eta_,ln_L_,ln_det_I_;
    fftw_complex    *FH_omp;

    using namespace std;
    //--- Some initialisation for computing BIC_c
    N = m - Ngaps; // Actual number of observations
    observations.get_values(m,&t_ori,&x_ori);
    observations.get_offsets(offsets);
    n_offsets = offsets.size();
    if (n_offsets<1) {
      cerr << "Not a single offset found! Cannot be good..." << endl;
      exit(EXIT_FAILURE);
    }
    ny = 2*m-1;  // array size needed for FFT operations
    nyc= ny/2+1; // arrayz size of complex values FFT arrays
    observations.get_t0t1(time0,time1); // start & end time of time series
    offset_index = designmatrix->get_offset_index(); // where are offsets in H
    column = offset_index + n_offsets - 1; // last offset column in H
    lambda = (time1-time0+1.0) / beta_spacing;
    cout << "lambda=" << lambda << endl;

    //--- For the OMP stuff I need Nthreads copy of the offsets
    off_omp.resize(Nthreads, vector<double>(n_offsets,0.0));
    for (i=0;i<Nthreads;i++) {
      for (j=0;j<n_offsets;j++) off_omp[i][j] = offsets[j];
    }

    //--- For the OMP stuff, I need to create various instances of several
    //    matrices and vectors to allow computations in parallel.
    try {
      H_omp          = new double[Nthreads*m*n];
      A1_omp         = new double[Nthreads*m*n];
      A2_omp         = new double[Nthreads*m*n];
      dummy_omp      = new double[Nthreads*ny];
      QA_omp         = new double[Nthreads*Ngaps*n];
      Qt_omp         = new double[Nthreads*Ngaps];
      C_theta_omp    = new double[Nthreads*n*n];
      C_thetaInv_omp = new double[Nthreads*n*n];
      z_omp          = new double[Nthreads*n];
      t1_omp         = new double[Nthreads*m];
      t2_omp         = new double[Nthreads*m];
      theta_omp      = new double[Nthreads*n];
      FH_omp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nyc*Nthreads);
    } 
    catch (bad_alloc()) {
      cerr << "AmmarGrag: need more memory!" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Fill all matrices and vectors
    for (i=0;i<Nthreads;i++) {
      cblas_dcopy(m,&H[m*column],1,&H_omp[i*m],1);
      cblas_dcopy(n,dummy,1,&z_omp[i*n],1);
      cblas_dcopy(m*n,A1,1,&A1_omp[i*m*n],1);
      cblas_dcopy(m*n,A2,1,&A2_omp[i*m*n],1);
      cblas_dcopy(Ngaps*n,QA,1,&QA_omp[i*Ngaps*n],1);
      cblas_dcopy(n*n,C_thetaInv,1,&C_thetaInv_omp[i*n*n],1);
    }

    //--- Create pivots (check which is larger n or Ngaps)
    k = n;
    if (Ngaps>n) k=Ngaps;
    ipiv = new int[k];
    for (i=0;i<k;i++) ipiv[i]=i+1;


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Loop over all possible offset locations, except first observation
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BIC_c[0] = LARGE;

    #pragma omp parallel for private(i,j,l,product,already_used,ln_det_I_,\
					     ln_L_,sigma_eta_,lnf_s,lnf_theta)
    for (i=1;i<m;i++) {

      //--- Which thread is this?
      l = omp_get_thread_num ();

      //--- Shift offset in design matrix H. Gaps are already  set to zero 
      //    but in this way I am sure that it always works correctly. Note
      //    that this subroutine assumes an offset is set at i=0.
      memset(&H_omp[l*m],0,i*sizeof(double));

      //--- Check if we already have an offset at this location or if we
      //    are dealing with a missing observations which needs to be skipped
      //    too.
      already_used = false;
      for (j=0;j<n_offsets-1;j++) {
        if (fabs(t[i]-off_omp[l][j])<EPS) {
          cout << "i=" << ", t=" << t[i] << endl;
          already_used=true;
        } 
      } 
      if (std::isnan(x_ori[i]) || already_used==true) {
        BIC_c[i] = LARGE;
      } else {

        //--- Set last offset to new location
        off_omp[l][n_offsets-1] = t[i];

        //--- Compute A1 & A2 for new design matrix H
        memset(&dummy_omp[l*ny+m],0,(m-1)*sizeof(double)); //--- zero padding 
        cblas_dcopy(m,&H_omp[l*m],1,&dummy_omp[l*ny],1);
        fftw_execute_dft_r2c(plan_forward,&dummy_omp[l*ny],&FH_omp[l*nyc]);

        #pragma omp critical
        {
          step2(1,&FH_omp[l*nyc],&A1_omp[l*m*n+m*column],
						    &A2_omp[l*m*n+m*column]);
        }

        //--- Update last column of inv(C_theta)
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,n,1,m,1.0,
			&A2_omp[l*m*n],m,&A2_omp[l*m*n+column*m],m, 
				   0.0,&C_thetaInv_omp[l*n*n+column*n],n);
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,n,1,m,1.0,
			&A1_omp[l*m*n],m,&A1_omp[l*m*n+column*m],m,
				  -1.0,&C_thetaInv_omp[l*n*n+column*n],n);
        z_omp[l*n+column]  = cblas_ddot(m,&A1_omp[l*m*n+column*m],1,y1,1);
        z_omp[l*n+column] -= cblas_ddot(m,&A2_omp[l*m*n+column*m],1,y2,1);

        //--- If required, add QA and Qy part to C_thetaInv and dummy      
        if (Ngaps>0) {
          //--- Update last column of QA
          cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		      Ngaps,1,m,1.0,G2,m,&A2_omp[l*m*n+column*m],m, 0.0,
				        &QA_omp[l*Ngaps*n+column*Ngaps],Ngaps);
          cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		      Ngaps,1,m,1.0,G1,m,&A1_omp[l*m*n+column*m],m,-1.0,
				        &QA_omp[l*Ngaps*n+column*Ngaps],Ngaps);
          dgetrs_(&Trans,&Ngaps,&one,M,&Ngaps,ipiv,
			        &QA_omp[l*Ngaps*n+column*Ngaps],&Ngaps,&info);
          //--- Adjust last column of C_thetaInv and z
          cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,n,1,Ngaps,-1.0,
	          &QA_omp[l*Ngaps*n],Ngaps,&QA_omp[l*Ngaps*n+column*Ngaps],
				  Ngaps,1.0,&C_thetaInv_omp[l*n*n+column*n],n);
          z_omp[l*n+column] -= cblas_ddot(Ngaps,
					&QA_omp[l*Ngaps*n+column*Ngaps],1,Qy,1);
        }

        //--- Compute C_theta
        cblas_dcopy(n*n,&C_thetaInv_omp[l*n*n],1,&C_theta_omp[l*n*n],1);
        dpotrf_(&Up,&n,&C_theta_omp[l*n*n],&n,&info);
			
        //--- Compute ln_det_I (The logarithm of the determinant of the
        //    Fisher Information matrix). det(C_theta^{-1}) = 1/det(C_theta)
        //    However, I = C_theta^{-1} so the inverses cancel.
        ln_det_I_ = 0.0;
        for (j=0;j<n;j++) {
          ln_det_I_ += 2.0*log(C_theta_omp[l*n*n + j+j*n]);
        }

        //--- Now inverse C_theta
        dpotri_(&Up,&n,&C_theta_omp[l*n*n],&n,&info);
      
        //--- Perform Least-Squares
        cblas_dsymv(CblasColMajor,CblasUpper,n,1.0,&C_theta_omp[l*n*n],n,
					&z_omp[l*n],1,0.0,&theta_omp[l*n],1);
								     

        //--- Compute residuals t1 and t2 = y1-A1*theta and y2-A2*theta
        cblas_dcopy(m,y1,1,&t1_omp[l*m],1);
        cblas_dcopy(m,y2,1,&t2_omp[l*m],1);
        cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,&A1_omp[l*m*n],m,
				        &theta_omp[l*n],1,-1.0,&t1_omp[l*m],1);
        cblas_dgemv(CblasColMajor,CblasNoTrans,m,n,1.0,&A2_omp[l*m*n],m,
					&theta_omp[l*n],1,-1.0,&t2_omp[l*m],1);

        product = cblas_ddot(m,&t1_omp[l*m],1,&t1_omp[l*m],1) - 
				cblas_ddot(m,&t2_omp[l*m],1,&t2_omp[l*m],1);
					
        if (Ngaps>0) {
          //--- Q_t = inv(M')*(G1'* t1 - G2'*t2) 
          cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,Ngaps,1,m,1.0,
			      G2,m,&t2_omp[l*m],m, 0.0,&Qt_omp[l*Ngaps],Ngaps);
          cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,Ngaps,1,m,1.0,
			      G1,m,&t1_omp[l*m],m,-1.0,&Qt_omp[l*Ngaps],Ngaps);
          dgetrs_(&Trans,&Ngaps,&one,M,&Ngaps,ipiv,
					       &Qt_omp[l*Ngaps],&Ngaps,&info);
          product -= cblas_ddot(Ngaps,&Qt_omp[l*Ngaps],1,&Qt_omp[l*Ngaps],1);
        }

        //--- Finally, compute sigma_eta
        sigma_eta_ = sqrt(product/static_cast<double>(m-Ngaps));

        //--- Fisher Information matrix was computed using unit covariance
        //    matrix. Update with right scaling.
        ln_det_I_ -= 2.0*n*log(sigma_eta_);

        //--- Compute Log-Likelihood
        ln_L_ = -0.5*(N*log(tpi) + ln_det_C + 2.0*N*log(sigma_eta_) + N);
        //cout << t[i] << "  " << ln_L << endl;

        lnf_s = 0.0;
        //--- Number of offsets follows Poisson distribution. I think I could
        //    move this section of code outside the loop over i to avoid 
        //    computing the same thing over and over but I've kept it here
        //    to keep the code readable.
        if (!std::isnan(beta_spacing)) {
          fact_k = 1.0;
          for (j=1;j<=n_offsets;j++) fact_k *= j;
          lnf_s = log(pow(lambda,n_offsets)*exp(-lambda)/fact_k);
        }

        //--- Size of offsets follows exponential distribution
        lnf_theta = 0.0;
        if (!std::isnan(beta_size)) {
          //--- prior for size of offsets     
          for (j=0;j<n_offsets;j++) {
            lnf_theta += -fabs(theta_omp[l*n+offset_index+j])/beta_size;
          }
          //cout << "lnf_theta= " << lnf_theta << endl;
        } 

        //--- Compute BIC_c (Nparam is defined in MLEBase)
        BIC_c[i] = (Nparam+n+1)*log(N) - 2.0*ln_L_ + ln_det_I_ 
		         - 2.0*lnf_s - 2.0*lnf_theta + extra_penalty*n_offsets;
        if (Kashyap==false) BIC_c[i] -= ln_det_I_;
      } 
    }

    //--- Free memory
    delete[] H_omp;
    delete[] A1_omp;
    delete[] A2_omp;
    delete[] dummy_omp;
    delete[] QA_omp;
    delete[] Qt_omp;
    delete[] C_theta_omp;
    delete[] C_thetaInv_omp;
    delete[] z_omp;
    delete[] t1_omp;
    delete[] t2_omp;
    delete[] theta_omp;
    delete[] ipiv;
    delete[] t_ori;
    delete[] x_ori;
    fftw_free (FH_omp);
  }
