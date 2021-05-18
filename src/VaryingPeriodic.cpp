/*! \filename VaryingPeriodic.cpp
 *  \author Machiel Bos
 *
 * Langbein used a bandpass noise model to deal with varying annual signals
 * He first defines the spectrum around the annual signal and from that
 * computes the impulse response coefficients. In Hector, the basis is
 * always the autocovariance function, not the spectrum or impulse response
 * coefficients. In 2018, we invented something similar using an AR(1)
 * amplitude variation, based on the work of Davis et al (2012).
 *
 * 4/5/2020  Santa Clara
 *
 * This script is part of Hector 1.9
 *
 * Hector is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * Hector is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hector. If not, see <http://www.gnu.org/licenses/>
 */
//==============================================================================
  #include "VaryingPeriodic.h"
  #include "Control.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstring>
  #include <cmath>

//==============================================================================
// Subroutines
//==============================================================================


//---!!-----------------------------------------------------------
  VaryingPeriodic::VaryingPeriodic(double T_0): pi(4.0*atan(1.0)),
                                                NaN(sqrt(-1.0))
//---!!-----------------------------------------------------------
  {
    Control   &control=Control::getInstance();

    using namespace std;
    //--- Remember the main frequency
    omega_0 = 2.0*pi/T_0;  // conert to rad/day
    try {
      control.get_string("PhysicalUnit",unit);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Check if phi is fixed or needs to be estimated
    try {
      phi_fixed = control.get_double("phi_varying_fixed");
      Nparam = 0;
    }
    catch (exception &e) {
      phi_fixed = NaN;
      Nparam = 1;
    }
  }



//-----------------------------------------------------------------
  void VaryingPeriodic::mactrick(int m, double *gamma_x, double *h)
//-----------------------------------------------------------------
  {
    double   sin_theta,cos_theta,*U,*V;
    int      k,k_old=0,k_new=1;

    using namespace std;
    //--- Define auxiliary arrays
    U = new double[2*m];
    V = new double[2*m];

    //--- define the generators u and v
    cblas_dcopy(m,gamma_x,1,U,1);
    cblas_dscal(m,1.0/sqrt(gamma_x[0]),U,1);
    cblas_dcopy(m,U,1,V,1);
    V[0] = 0.0;
    h[m-1] = U[m-1];

    for (k=1;k<m;k++) {
      if (U[(k-1) + m*k_old]<1.0e-4) {
        cout << "mactrick: small U[k-1 + m*k_old]=" << U[k-1 + m*k_old] << endl;
        exit(EXIT_FAILURE);
      }
      sin_theta = V[k + m*k_old]/U[(k-1) + m*k_old]; // Eq. (3.3a)
      cos_theta = sqrt(1.0-sin_theta*sin_theta);     // Eq. (3.3b)
      if (cos_theta<1.0e-4) {
        cout << "mactrick: cos_theta=" << cos_theta << endl;
        exit(EXIT_FAILURE);
      }

      //--- Eq. (3.7a)
      cblas_dcopy(m-k,&V[k + m*k_old],1,&V[k + m*k_new],1);
      cblas_daxpy(m-k,-sin_theta,&U[k-1 + m*k_old],1,&V[k + m*k_new],1);
      cblas_dscal(m-k,1.0/cos_theta,&V[k + m*k_new],1);

      //-- Eq. (3.7b)
      cblas_dcopy(m-k,&U[(k-1) + m*k_old],1,&U[k + m*k_new],1);
      cblas_dscal(m-k,cos_theta,&U[k + m*k_new],1);
      cblas_daxpy(m-k,-sin_theta,&V[k + m*k_new],1,&U[k + m*k_new],1);
      h[m-1-k] = U[m-1 + m*k_new];

      k_old = k_new;
      k_new = 1 - k_new;
    }

    //--- Free memory
    delete[] U;
    delete[] V;
  }



/*! The covariance matrix is a unit matrix
 */
//---------------------------------------------------------------------------
  void VaryingPeriodic::get_covariance(double *param, int m, double *gamma_x)
//---------------------------------------------------------------------------
  {
    int       i;
    double    frac,phi;

    using namespace std;
    if (std::isnan(phi_fixed)==true) {
      phi = param[0];
    } else {
      phi = phi_fixed;
    }

    frac = 0.5/(1 - phi*phi);
    gamma_x[0] = frac;
    for (i=1;i<m;i++) {
      frac *= phi;
      gamma_x[i] = frac*cos(omega_0 * static_cast<double>(i));
    }
  }



/*! Nothing to show
 */
//----------------------------------------------------------------------
  void VaryingPeriodic::show(double *param, double *error, double sigma)
//----------------------------------------------------------------------
  {
    JSON  &json = JSON::getInstance();

    using namespace std;
    cout << "sigma     = " << sigma << " " << unit << endl;
    json.write_double("sigma",sigma);

    if (Nparam==1) {
      cout << "phi       = " << param[0] << " +/- " << error[0] << endl;
      json.write_double("phi",param[0]);
    } else {
      cout << "phi       = " << phi_fixed << " (fixed)" << endl;
      json.write_double("phi",phi_fixed);
    }
  }



/*! There is no penalty if phi_fixed is set
 */
//------------------------------------------------------
  double VaryingPeriodic::compute_penalty(double *param)
//------------------------------------------------------
  {
    const double LARGE=1.0e5;
    double       phi,penalty=0.0;

    if (Nparam==1) {
      phi = param[0];
      if (phi<0.5) {
        param[0] = 0.5;
        penalty = (0.5-phi)*LARGE;
      } else if (phi>0.999999) {
        param[0] = 0.999999;
        penalty = (phi-0.999999)*LARGE;
      }
    } else {
      penalty = 0.0;
    }

    return penalty;
  }



/*! Nothing can be set or needs to be stored for white noise
 */
//----------------------------------------------------------------
  void VaryingPeriodic::set_noise_parameters(double *params_fixed)
//----------------------------------------------------------------
  {
    using namespace std;
    if (Nparam==1) {
      cout << "Enter value for phi: ";
      cin >> phi_fixed;
    }
  }



/*! Compute PSD for given frequency
 */
//------------------------------------------------
  double VaryingPeriodic::compute_G(double lambda)
//------------------------------------------------
  {
    double   phi;

    using namespace std;
    if (std::isnan(phi_fixed)==true) {
      cerr << "phi_fixed is not set!!! ERROR" << endl;
      exit(EXIT_FAILURE);
    }
    phi = phi_fixed;

    return 1.0/2.0 * (1.0/(1.0 - 2.0*phi*cos(lambda+omega_0) + phi*phi) +
                      1.0/(1.0 - 2.0*phi*cos(lambda-omega_0) + phi*phi));
  }



/*! Compute impulse response: h
 */
//----------------------------------------------------------------
  void VaryingPeriodic::compute_impulse_response(int m, double* h)
//----------------------------------------------------------------
  {
    int    i;
    double I,params_fixed[1],*gamma_x;

    //--- For this particular noise model, I need gamma_x and perform Cholesky
    gamma_x = new double[m];

    if (Nparam==1) {
      set_noise_parameters(params_fixed);
    }
    //--- Eq. (72)
    get_covariance(params_fixed,m,gamma_x);
    mactrick(m,gamma_x,h);

    delete[] gamma_x;
  }
