/*! \file   MLEBase.cpp
 *  \author Machiel Bos
 *
 * The class that prepares the arrays before MLE can be computed using
 * AmmarGrag or FullCov.
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
 *  along with Hector. If not, see <http://www.gnu.org/licenses/>
 *
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
  MLEBase::MLEBase(void) : tpi(8.0*atan(1.0)),
                           NaN(sqrt(-1.0))
//---!!---------------------------------------
  {
    Control        &control = Control::getInstance();
    Observations   &observations=Observations::getInstance();
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    NoiseModel     &noisemodel=NoiseModel::getInstance(); 

    using namespace std;
    //--- Prior probability of size of offset
    try {
      beta_size  = control.get_double("beta_size");
      if (beta_size<0.0) {
        cerr << "Only positive beta's are allowed! : " << beta_size << endl;
        exit(EXIT_FAILURE);
      }
      cout << "beta_size: " << beta_size << endl;
    }
    catch (exception &e) {
      beta_size    = NaN;
    }

    //--- Prior probability of spacing of offset
    try {
      beta_spacing = 365.25*control.get_double("beta_spacing");
      if (beta_spacing<0.0) {
        cerr << "Only positive beta's are allowed! : " << beta_spacing << endl;
        exit(EXIT_FAILURE);
      }
      cout << "beta_spacing: " << beta_spacing/365.25 << endl;
    }
    catch (exception &e) {
      beta_spacing = NaN;
    }

    //--- Use Fisher Information matrix or not in BIC_c
    try {
      Kashyap = control.get_bool("Kashyap");
      cout << "use Kashyap: " << Kashyap << endl;
    }
    catch (exception &e) {
      Kashyap = true;
    }

    //--- What is the extra penalty scale factor for each added parameter?
    try {
      extra_penalty = control.get_double("BIC_c_ExtraPenalty");
      cout << "Penalty in BIC_c for each parameter: " << extra_penalty << endl;
    }
    catch (exception &e) {
      extra_penalty = 2.0;
    }

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
    theta      = new double[n];
    C_theta    = new double[n*n];
    C_thetaInv = new double[n*n];
    memset(theta,0,(n)*sizeof(double));     // put theta to zero
    memset(C_theta,0,(n*n)*sizeof(double)); // put C_theta to zero
  }



/*! Free up Memory 
 */
//-----------------------
  MLEBase::~MLEBase(void)
//-----------------------
  {
    if (t!=NULL)          delete[]  t;
    if (x!=NULL)          delete[]  x;
    if (H!=NULL)          delete[]  H;
    if (F!=NULL)          delete[]  F;
    if (theta!=NULL)      delete[]  theta;
    if (C_theta!=NULL)    delete[]  C_theta;
    if (C_thetaInv!=NULL) delete[]  C_thetaInv;
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
    JSON           &json=JSON::getInstance();

    using namespace std; 
    //---- Already show value of sigma_eta
    cout << "STD of the driving noise: " << sigma_eta << endl;
    json.write_double("driving_noise",sigma_eta);
 
    error = new double[n];
    for (i=0;i<n;i++) error[i] = sigma_eta*sqrt(C_theta[i+n*i]);
    designmatrix->show_results(theta,error);
    designmatrix->compute_xhat(theta);
    observations.save_mom(true);

    delete[] error;
  }



/*! Computed the likelihood function. 
 */
//---------------------------------------
  double MLEBase::compute(double *param_)
//---------------------------------------
  {
    int           N;
    double        l,penalty,*param=NULL;
    NoiseModel    &noisemodel=NoiseModel::getInstance();

    using namespace std;
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
    compute_LeastSquares(param);

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
    vector<double>  offsets;
    int             i,j,k,n_offsets,offset_index,N;
    double          t0,t1,lambda,fact_k;
    double          lnf_s,lnf_theta;
    DesignMatrix    *designmatrix = DesignMatrix::getInstance();
    Observations    &observations = Observations::getInstance();

    //--- How many parameters are estimated? 
    k = Nparam + n + 1;

    //--- Compute probability of the size and spacing of offsets
    observations.get_offsets(offsets);
    n_offsets = offsets.size();
    lnf_s     = 0.0;
    //--- Number of offsets follows poisson distribution, also need to be
    //    computed for zero offsets!
    if (!std::isnan(beta_size)) {
      observations.get_t0t1(t0,t1);
      offset_index = designmatrix->get_offset_index();

      lambda = (t1-t0+1.0) / beta_spacing;
      fact_k = 1.0;
      for (j=1;j<=n_offsets;j++) fact_k *= j;
      lnf_s = log(pow(lambda,n_offsets)*exp(-lambda)/fact_k);

    }

    lnf_theta = 0.0;
    //--- prior for size of offsets     
    if (n_offsets>0 && !std::isnan(beta_spacing)) {
      for (j=0;j<n_offsets;j++) {
        lnf_theta += -fabs(theta[offset_index+j])/beta_size;
      }
    }

    //--- Compute log(likelihood) at minimum
    N     = m - Ngaps;
    ln_L  = -compute(param);
    AIC   = 2.0*k - 2.0*ln_L;
    BIC   = k*log(N) - 2.0*ln_L;
    BIC_tp= k*log(N/tpi) - 2.0*ln_L;
    BIC_c = k*log(N) - 2.0*ln_L + ln_det_I - 2.0*lnf_s 
		                   - 2.0*lnf_theta + extra_penalty*n_offsets;
    if (Kashyap==false) BIC_c -= ln_det_I;
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
      //for (j=0;j<n;j++) printf("%7.3lf  ",A[i+m*j]);
      //cout << endl;
      for (j=0;j<n;j++) fp << setprecision(5) << setw(10) << A[i+m*j] << "  ";
      fp << endl;
    }
    fp.close();
  }
