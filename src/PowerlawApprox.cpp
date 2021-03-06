/*! \file    PowerlawApprox.cpp
 *  \author  Machiel Bos
 *
 * The most common noise model for GPS observations is power-law + white
 * noise. I switched to using d (=alpha/2) and phi. If d<0.5, then there
 * exist closed formula's for the covariance matrix. However, for d>=0.5
 * the variance is infinite. However, the growth of the variance is slow.
 * The trick I use is to let the noise start to develop at some time
 * before the first observation and let the noise variance grow up to the
 * last observed epoch. The column and row of the covariance matrix for the
 * last epoch is used to create a Toeplitz matrix, see Bos et al. (2013).
 *
 *  This script is part of Hector 1.7.2
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
  #include "PowerlawApprox.h"
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <iomanip>
  #include <cmath>
  #include <cstdlib>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!--------------------------------------------------------------
  PowerlawApprox::PowerlawApprox(double d_fixed_) : pi(4.0*atan(1.0))
//---!!--------------------------------------------------------------
  {
    Control       &control=Control::getInstance();

    using namespace std;
    //--- Check if we need to estimate the spectral index or not
    if (!std::isnan(d_fixed_)) {
      estimate_spectral_index = false;
      d_fixed = d_fixed_;
    } else {
      estimate_spectral_index = true;
    }

    //--- How many days before t0 did the noise start to develop?
    try {
      ms = control.get_int("TimeNoiseStart"); // number of days before t0
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Return the covariance matrix (only first column)
 *
 *  \param[in]   param   : param[0]  = d,  param[1]  = phi
 *  \param[in]   m       : length of gamma_x
 *  \param[out]  gamma_x : covariance matrix (first column)  
 */
//--------------------------------------------------------------------------
  void PowerlawApprox::get_covariance(double *param, int m, double *gamma_x) 
//--------------------------------------------------------------------------
  {
    int                 i,ny;
    const double        TINY=1.0e-6;
    double              Scale,I,d,d_max;
    double              *C=NULL,*h1=NULL,*h2=NULL;
    fftw_plan           plan_forward,plan_backward;
    fftw_complex        *F_h1=NULL,*F_h2=NULL,*F_C=NULL;

    using namespace std;
    //--- An approximation is an approximation. Set d_max
    d_max = 0.85;

    //--- spectral index is stored in first element
    if (estimate_spectral_index==true) {
      d = param[0];
    } else {
      d = d_fixed;
    }

    //--- Sanity check
    if (d-TINY>=d_max || d+TINY<=-0.5) {
      cerr << "PowerlawApprox: d is outside its range: " << d << endl;
      exit(EXIT_FAILURE);
    }

    ny   = 2*(ms+m)-1;  // We have 2*(ms+m)-1 values for circulant matrix
    h1   = new double[ny]; //--- Normal array with h coeffients
    h2   = new double[ny]; //    same as h1 but in reverse order
    C    = new double[ny]; //    just an array but represents last column of C
    F_h1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (ny/2+1));
    F_h2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (ny/2+1));
    F_C  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (ny/2+1));

    //--- Remember that we have have power-law + white.
    h1[ms+m-1] = 1.0;

    //--- Apply Kasdin recursion formula (fill 1 from m to ms+m
    for (i=1,I=1.0;i<ms+m;i++,I+=1.0) {
      h1[ms+m-1-i] = (d + I - 1.0)*h1[ms+m-i]/I;
    }

    //--- Reverse h1 into h2
    cblas_dcopy(ms+m,h1,1,&h2[0],-1);

    //--- pad rest (from m+ms until (m+ms)+(m+ms-1) with zeros to make circulant
    for (i=ms+m;i<2*(ms+m)-1;i++) {
      h1[i]=h2[i]=0.0;
    }
#ifdef DEBUG
    cout << "ms=" << ms << endl;
    for (i=0;i<(ms+m);i++) {
      cout << "i=" << i << ",  h1=" << h1[i] << ", h2=" << h2[i] << endl;
    }
#endif

    //--- Create FFT forward & backward plan
    plan_forward  = fftw_plan_dft_r2c_1d(ny,   h1, F_h1,FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(ny, F_h1,   h1,FFTW_ESTIMATE);
    Scale  = 1.0/static_cast<double>(ny);
 
    //--- Compute forwards FFT for h1 and h2
    fftw_execute_dft_r2c(plan_forward,h1,F_h1);
    fftw_execute_dft_r2c(plan_forward,h2,F_h2);

    //--- Perform C = h1*h2 termwise 
    for (i=0;i<(ny/2+1);i++) {
      F_C[i][0] = F_h1[i][0]*F_h2[i][0] - F_h1[i][1]*F_h2[i][1];
      F_C[i][1] = F_h1[i][0]*F_h2[i][1] + F_h1[i][1]*F_h2[i][0];
    }

    //--- Compute backward FFT
#ifdef DEBUG
    cout << "F_C before iFFT=" << endl;
    for (i=0;i<(ny/2+1);i++) {
      cout << "i=" << i << ",  C=" << F_C[i][0] << ", " << F_C[i][1] << endl;
    }
#endif
    fftw_execute_dft_c2r(plan_backward,F_C,C);
#ifdef DEBUG
    cout << "C after iFFT=" << endl;
    for (i=0;i<ny;i++) {
      cout << "i=" << i << ",  C=" << C[i] << endl;
    }
#endif

    //--- Form covariance matrix by copying power-law covariance (last column
    //    of C, into gamma_x. The larger ms, the better this column represents
    //    the toeplitz nature of C
    cblas_dscal(m,Scale,&C[ms],1);
    cblas_dcopy(m,&C[ms],-1,gamma_x,1);

    //--- Free memory
    delete[]  h1;
    delete[]  h2;
    delete[]  C;

    //--- free FFTW3 variables 
    fftw_destroy_plan (plan_forward);
    fftw_destroy_plan (plan_backward);
    fftw_free ( F_h1 );
    fftw_free ( F_h2 );
    fftw_free ( F_C  );
  }



/*! Show noise parameters. 
 */
//---------------------------------------------------------------------
  void PowerlawApprox::show(double *param, double *error, double sigma)
//---------------------------------------------------------------------
  {
    double        d,T;
    Observations  &observations=Observations::getInstance();

    using namespace std;
    T = 1.0/(365.25*24.0*3600.0*observations.get_fs()); // T in yr

    if (estimate_spectral_index==true) {
      d = param[0];
    } else {
      d = d_fixed;
    }

    cout << "sigma     = " << sigma/pow(T,0.5*d)
                          << " " << unit << "/yr^" << 0.5*d << endl;
    cout << fixed << setprecision(4); 
    cout << "d         = " << d << " +/- " << error[0] << endl;
    cout << "kappa     = " << -2*d << " +/- " << 2*error[0] << endl;
  }



/*! Number of parameters to estimate is 1: alpha = 2*d
 */
//-------------------------------------------
  int PowerlawApprox::get_Nparam(void)
//-------------------------------------------
  {
    if (estimate_spectral_index==true) {
      return 1;
    } else {
      return 0;
    }
  }



/*! alpha has a valid range between -1 and 1 -> d element of [-0.5:0.5]
 */
//-----------------------------------------------------
  double PowerlawApprox::compute_penalty(double *param)
//-----------------------------------------------------
  {
    double   penalty=0.0,LARGE=1.0e8;

    if (estimate_spectral_index==true) {
      //--- d
      if (param[0]>0.85) {
        penalty += (param[0]-0.85)*LARGE;
        param[0] = 0.85;
      } 
      if (param[0]<-0.5) {
        penalty += (-0.5-param[0])*LARGE;
        param[0] = -0.5;
      }
    }

    return penalty;
  }


/*! Not implemented
 */
//---------------------------------------------------------------
  void PowerlawApprox::set_noise_parameters(double * /*params_fixed*/)
//---------------------------------------------------------------
  {
    using namespace std;
    cerr << "Please use pure Powerlaw class!" << endl;
    exit(EXIT_FAILURE);
  }



/*! Compute PSD for given frequency
 */
//-----------------------------------------------
  double PowerlawApprox::compute_G(double /*lambda*/)
//-----------------------------------------------
  {
    using namespace std;
    cerr << "Please use pure Powerlaw class!" << endl;
    exit(EXIT_FAILURE);
  }



/*! Compute impulse response: h
 */
//---------------------------------------------------------------
  void PowerlawApprox::compute_impulse_response(int /*m*/, double* /*h*/)
//---------------------------------------------------------------
  {
    using namespace std;
    cerr << "Please use pure Powerlaw class!" << endl;
    exit(EXIT_FAILURE);
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
