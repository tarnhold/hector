/*! \file   Matern.cpp
 *  \author Machiel Bos
 *
 * This class implements the Matérn noise model. Eqs. are taken from
 * Lilly et al. (2017) "Fractional Brownian motion, the Matérn process, and 
 * stochastic modeling of turbulent dispersion", Nonlinear Processes in 
 * Geophysics, 24: 481-514 
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
  #include "Matern.h"
  #include <iostream>
  #include <ostream>
  #include <cmath>
  #include <cstdlib>
  #include <cstdio>
  #include <cstring>

  #define TINY 1.0e-12
//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!---------------------------------------------
  Matern::Matern(double d_fixed) : NaN(sqrt(-1.0)),
                                   pi(4.0*atan(1.0))
//---!!---------------------------------------------
  {
    Control  &control=Control::getInstance();

    using namespace std;
    //--- As always, check if Physical Unit is provided
    try {
      control.get_string("PhysicalUnit",unit);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    //--- Check if spectral index is fixed or needs to be estimated
    if (std::isnan(d_fixed)==true) {
      alpha_fixed = NaN;
      Nparam = 1;
    } else {
      alpha_fixed = 2.0*d_fixed;
      cout << "alpha fixed to : " << alpha_fixed << endl;
      Nparam = 0;
    }
    
    //--- Check if lambda is fixed or needs to be estimated
    try {
      lambda_fixed = control.get_double("lambda_fixed");
    }
    catch (exception &e) {
      lambda_fixed = NaN;
      Nparam += 1;
    }
  }



//--------------------------------------------------------
  void Matern::mactrick(int m, double *gamma_x, double *h)
//--------------------------------------------------------
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



/*! Compute covariance vector
 */
//------------------------------------------------------------------
  void Matern::get_covariance(double *param, int m, double *gamma_x)
//------------------------------------------------------------------
  {
    using namespace std;
    const double  threshold=1.0e3,EPS=1.0e-6;
    int           i;
    double        tau,c0,K_nu,K_nu_approx,d,alpha,lambda;

    //--- extract parameters to readable variables
    if (Nparam==0) {
      alpha  = alpha_fixed;
      lambda = lambda_fixed;
    } else if (Nparam==1) {
      if (std::isnan(alpha_fixed)==true) {
        d      = param[0];
        lambda = lambda_fixed;
      } else {
        d      = 0.5*alpha_fixed;
        lambda = param[0];
      }
    } else {
      d      = param[0];
      lambda = param[1];
    }
    alpha = 2.0*d;

    //--- Constant
    c0 = 2.0/(tgamma(alpha-0.5) * pow(2.0, alpha-0.5));

    //--- Use Modified Bessel Function, second order
    i = 0;
    tau = static_cast<double>(i);
    while (i<m && lambda*tau<threshold) {
      if (fabs(lambda*tau)<EPS) {
        //--- Eq. (62)
        gamma_x[i] = 1.0 - pow(lambda*tau/2.0,2.0*alpha-1.0)*tgamma(1.5-alpha)/
						          tgamma(alpha+0.5);
      } else {
        //--- Eq. (60)
        gamma_x[i] = c0 * pow(lambda*tau,alpha-0.5) * \
			boost::math::cyl_bessel_k(alpha-0.5,lambda*tau);
      }
      i++;
      tau = static_cast<double>(i);
    }

    if (i<m) {
      //--- Check approximation
      tau = static_cast<double>(i);
      K_nu = boost::math::cyl_bessel_k(alpha-0.5,lambda*tau);
      K_nu_approx = sqrt(0.5*pi) * pow(lambda*tau,-0.5) *exp(-lambda*tau);

      if (fabs((K_nu - K_nu_approx)/K_nu)>EPS) {
        cerr << "Oops" << endl;
        cerr << "K_nu=" << K_nu << endl;
        cerr << "K_nu_approx=" << K_nu_approx << endl;
        exit(EXIT_FAILURE);
      }
    }

    //--- Use approximation for the rest of the autocovariance function
    while (i<m) {
      tau = static_cast<double>(i);
      //--- Eq. (61)
      K_nu_approx = sqrt(0.5*pi) * exp(-lambda*tau);
      gamma_x[i] = c0 * pow(lambda*tau,alpha-1.0) * K_nu_approx;
      i++;
    }
    //cout << "alpha=" << alpha << ", lambda=" << lambda << endl;
  }



/*! Show phi and fractional difference values
 */
//-------------------------------------------------------------
  void Matern::show(double *param, double *error, double sigma)
//-------------------------------------------------------------
  {
    double        d,T;
    Observations  &observations=Observations::getInstance();
    JSON          &json = JSON::getInstance();

    using namespace std;
    T = 1.0/(365.25*24.0*3600.0*observations.get_fs()); // T in yr

    if (std::isnan(alpha_fixed)) d = param[0];
    else                         d = 0.5*alpha_fixed;
    
    cout << "sigma     = " << sigma/pow(T,0.5*d)
             			<< " " << unit << "/yr^" << 0.5*d << endl;   
     
    if (Nparam==0) {
      printf("d         = %8.4lf (fixed)\n",0.5*alpha_fixed);
      printf("kappa     = %8.4lf (fixed)\n",-1.0*alpha_fixed);
      printf("lambda     = %le (fixed)\n",lambda_fixed);
      json.write_double("d",0.5*alpha_fixed);
      json.write_double("kappa",-1.0*alpha_fixed);
      json.write_double("lambda",lambda_fixed);
    } else if (Nparam==1) {
      if (std::isnan(alpha_fixed)==true) {
        printf("d         = %8.4lf +/- %6.4lf\n",param[0],error[0]);
        printf("kappa     = %8.4lf +/- %6.4lf\n",-2.0*param[0],2.0*error[0]);
        printf("lambda    = %le (fixed)\n",lambda_fixed);
        json.write_double("d",param[0]);
        json.write_double("kappa",-2.0*param[0]);
        json.write_double("lambda",lambda_fixed);
      } else {
        printf("d         = %8.4lf (fixed)\n",0.5*alpha_fixed);
        printf("kappa     = %8.4lf (fixed)\n",-1.0*alpha_fixed);
        printf("lambda    = %le +/- %le (fixed)\n",param[0],error[0]);
        json.write_double("d",0.5*alpha_fixed);
        json.write_double("kappa",-1.0*alpha_fixed);
        json.write_double("lambda",param[0]);
      }
    } else {
      printf("d         = %8.4lf +/- %6.4lf\n",param[0],error[0]);
      printf("kappa     = %8.4lf +/- %6.4lf\n",-2.0*param[0],2.0*error[0]);
      printf("lambda    = %le +/- %le\n",param[1],error[1]);
      json.write_double("d",param[0]);
      json.write_double("kappa",-2.0*param[0]);
      json.write_double("lambda",param[1]);
    }
  }



/*! Compute penalty function is parameters get out of bound
 */
//---------------------------------------------
  double Matern::compute_penalty(double *param)
//---------------------------------------------
  {
    double        penalty=0.0;
    const double  LARGE=1.0e5;

    using namespace std;
    if (Nparam==0) {
      penalty = 0.0;
    } else {
      if (Nparam==2 or (Nparam==1 and std::isnan(lambda_fixed)==false)) {
        //--- alpha is between 0.5 and 1.5
        if (param[0]>0.7499) {
          penalty += (param[0]-0.7499)*LARGE;
          param[0] = 0.7499;
        } else if (param[0]<0.251) {
          penalty += (0.251-param[0])*LARGE;
          param[0] = 0.251;
        }
      }
      if (Nparam==2) {
        //--- lambda>0  this is second parameter to estimate
        if (param[1]<1.0e-6) {
          penalty += abs(1.0e-6-param[1])*1.0e10;
          param[1] = 1.0e-6;
        } 
      } 
      if (Nparam==1 and std::isnan(alpha_fixed)==false) {
        if (param[0]<1.0e-8) {
          penalty += abs(1.0e-8-param[0])*1.0e12;
          param[0] = 1.0e-8;
        } 
      } 
    }

    return penalty;
  }



/*! Set noise parameters for PSD computation.
 *
 * If phi_fixed is set, only d can vary and Nparam==1. 
 * If d_fixed is also set, nothing can vary and Nparam==0.
 * Otherwise, set both d and phi
 */
//-------------------------------------------------------
  void Matern::set_noise_parameters(double *params_fixed)
//-------------------------------------------------------
  {
    double    d;
    using namespace std;
    if (Nparam>=1) {
      if (std::isnan(alpha_fixed)==true) {
        cout << "Enter parameter value of d: ";
        cin >> d;
        alpha_fixed = 2.0*d;
      }
    }
    if (Nparam==2) {
      if (std::isnan(lambda_fixed)==true) {
        cout << "Enter parameter value of lambda: ";
        cin >> lambda_fixed;
      }
    }
  }



/*! Compute PSD for given frequency
 */
//--------------------------------------
  double Matern::compute_G(double omega)
//--------------------------------------
  {
    double   c_alpha;

    using namespace std;
    //--- Eq. (56)
    c_alpha = tgamma(0.5)*tgamma(alpha_fixed-0.5)/(2.0*pi*tgamma(alpha_fixed));
  
    /* 
    cout << "c_alpha=" << c_alpha << endl;
    cout << "alpha_fixed=" << alpha_fixed << endl;
    cout << "lambda_fixed=" << lambda_fixed << endl;
    cout << "omega=" << omega << endl;
    */
 
    //--- Eq. (59)
    return pow(lambda_fixed,2.0*alpha_fixed-1.0)/c_alpha * 
	     	1.0/pow(pow(omega,2.0) + pow(lambda_fixed,2.0),alpha_fixed);
  }



/*! Compute impulse response: h
 */
//-------------------------------------------------------
  void Matern::compute_impulse_response(int m, double* h)
//-------------------------------------------------------
  {
    int    i;
    double I,*params_fixed=NULL,*gamma_x;

    //--- For this particular noise model, I need gamma_x and perform Cholesky
    gamma_x = new double[m];

    //--- get values of noise parameters
    if (Nparam>0) params_fixed = new double[Nparam];
    set_noise_parameters(params_fixed);
    params_fixed[0] = 0.5*alpha_fixed;
    params_fixed[1] = lambda_fixed;

    //--- Eq. (72)
    get_covariance(params_fixed,m,gamma_x);
    mactrick(m,gamma_x,h);

    if (params_fixed!=NULL) delete[] params_fixed;
    delete[] gamma_x;
  }
