/*! \file   Powerlaw.cpp
 *  \author Machiel Bos
 *
 * The most common noise model for GPS observations is power-law + white
 * noise. This class provides the power-law model. I switched to 
 * using d (=alpha/2) and phi.
 *
 * 7/10/2012  Santa Clara
 */
//==============================================================================
  #include "Powerlaw.h"
  #include <iostream>
  #include <iomanip>
  #include <ostream>
  #include <cmath>
  #include <cstdlib>

//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!--------------------------------------------------
  Powerlaw::Powerlaw(double d_fixed_) : pi(4.0*atan(1.0))
//---!!--------------------------------------------------
  {
    if (isnan(d_fixed_)==false) {
      estimate_spectral_index = false;
      d_fixed = d_fixed_;
    } else {
      estimate_spectral_index = true;
    }
  }



/*! Return the covariance matrix (only first column)
 *
 *  \param[in]   param   : param[0]  = d,  param[1]  = phi
 *  \param[in]   m       : length of gamma_x
 *  \param[out]  gamma_x : covariance matrix (first column)  
 */
//--------------------------------------------------------------------
  void Powerlaw::get_covariance(double *param, int m, double *gamma_x)
//--------------------------------------------------------------------
  {
    int            j;
    const double   TINY=1.0e-6;
    double         d,d_max,alpha,tau;
    Control        *control=Control::getInstance();

    using namespace std;
    //--- spectral index is stored in first element
    if (estimate_spectral_index==true) {
      d = param[0];
    } else {
      d = d_fixed;
    }
    
    if (control->get_bool("firstdifference")==true) { 
      d_max = 1.499;
      alpha = 2*d-2;
    } else {
      d_max = 0.499;
      alpha = 2*d;
    }

    //--- Sanity check
    if (d-TINY>=d_max || d+TINY<=-0.5) {
      cerr << "d is outside its possible range: " << d << endl;
      exit(EXIT_FAILURE);
    }
 
    //--- Form covariance matrix
    gamma_x[0] = exp(gamma(1.0-alpha))/pow(exp(gamma(1.0-0.5*alpha)),2.0);
					
    for (j=0;j<m-1;j++) {
      tau = static_cast<double>(j+1);
      gamma_x[j+1] = (0.5*alpha+tau-1.0)*gamma_x[j]/(tau-alpha/2.0);
    }
  }



/*! Show noise parameters. Since I use Simon's trick of factoring out the
 *  innovation noise variance out of the covariance matrix, only the angle
 *  phi and spectral index alpha are known. I show cos(phi) and sin(phi)
 *  to facilitate computing sigma_pl and sigma_w.
 */
//-------------------------------------------------
  void Powerlaw::show(double *param, double *error)
//-------------------------------------------------
  {
    double   d,phi;

    using namespace std;
    if (estimate_spectral_index==true) {
      d = param[0];
    } else {
      d = d_fixed;
    }

    cout << fixed << setprecision(4);
    cout << "d = " << d << " +/- " << error[0] << endl;
  }



/*! Number of parameters to estimate is 1: d
 */
//------------------------------
  int Powerlaw::get_Nparam(void)
//------------------------------
  {
    if (estimate_spectral_index==true) {
      return 1;
    } else {
      return 0;
    }
  }



/*! alpha has a valid range between -1 and 1 -> d element of [-0.5:0.5]
 */
//-----------------------------------------------
  double Powerlaw::compute_penalty(double *param)
//-----------------------------------------------
  {
    double   penalty=0.0,LARGE=1.0e8;
    Control  *control=Control::getInstance();

    using namespace std;
    if (estimate_spectral_index==true) {
      //--- d
      if (control->get_bool("firstdifference")==true) {
        if (param[0]>1.499) {
          penalty += (param[0]-1.499)*LARGE;
          param[0] = 1.499;
        }
      } else {
        if (param[0]>0.499) {
          penalty += (param[0]-0.499)*LARGE;
          param[0] = 0.499;
        }
      } 
      if (param[0]<-0.5) {
        penalty += (-0.5-param[0])*LARGE;
        param[0] = -0.5;
      }
#ifdef DEBUG
      cout << "Powerlaw: out=" << param[0] << ", penalty=" << penalty << endl;
#endif
   }

   return penalty;
  }



/*! Store angle and spectral index
 */
//------------------------------
  void Powerlaw::setup_PSD(void)
//------------------------------
  {
    using namespace std;
    if (estimate_spectral_index==true) {
      cout << "Enter value of fractional difference d:";
      cin >> param_PSD[0];
    } else {
      param_PSD[0] = d_fixed;
    }
  }



/*! Compute PSD for given frequency
 */
//-----------------------------------------
  double Powerlaw::compute_G(double lambda)
//-----------------------------------------
  {
    double    d;

    d = param_PSD[0];

    return 1.0/pow(2.0*sin(0.5*lambda),2.0*d);
  }

