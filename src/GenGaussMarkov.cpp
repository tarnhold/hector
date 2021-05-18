/*! \file   GenGaussMarkov.cpp
 *  \author Machiel Bos
 *
 * This class implements the Generalized Gauss-Markov noise model that is
 * preferred by John Langbein and Simon Williams.
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
  #include "GenGaussMarkov.h"
  #include <boost/math/special_functions/next.hpp>
  #include <boost/math/special_functions/hypergeometric_pFq.hpp>
  #include <boost/multiprecision/cpp_int.hpp> 
  #include <boost/multiprecision/cpp_bin_float.hpp> 
  #include <iostream>
  #include <ostream>

  #define TINY 1.0e-15
  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!------------------------------------------
  GenGaussMarkov::GenGaussMarkov(double d_fixed_) : NaN(sqrt(-1.0))
//---!!------------------------------------------
  {
    Control  &control=Control::getInstance();

    using namespace std;
    d_fixed = d_fixed_;
    try {
      control.get_string("PhysicalUnit",unit);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    if (std::isnan(d_fixed)==false) {
      //--- If we set d fixed, then phi is automatically fixed too. The idea
      //    is to have a small value so that the power-law does not flatten 
      //    out too soon.
      try {
        //--- phi is actually stored as 1-phi (confusing...)
        phi_fixed = control.get_double("GGM_1mphi");
      }
      catch (exception &e) {
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
      }

      //--- With everything fixed, the number of parameters to estimate is 0
      Nparam = 0;
    } else {
      try {
        //--- phi is actually stored as 1-phi (confusing...)
        phi_fixed = control.get_double("GGM_1mphi");
        Nparam = 1;
      }
      catch (exception &e) {
        phi_fixed = NaN;
        Nparam = 2;
      }
    }
  }

 

/*! Compute backward recursion
 */
//--------------------------------------------------------------------------
  double GenGaussMarkov::backward(double a, double b, double c, 
                                             double z, double F, double Fp1)
//--------------------------------------------------------------------------
  {
    return ((1.0-c+(b-a)*z)*F + (a*(c-b)*z)*Fp1/c)/(1.0-c);
  }



/*! Compute covariance vector
 */
//--------------------------------------------------------------------------
  void GenGaussMarkov::get_covariance(double *param, int m, double *gamma_x)
//--------------------------------------------------------------------------
  {
    using namespace std;
    int              i,k;
    double           scale,tau,phi,a,b,c,d,z,*_2F1=NULL,Fm1,F,Fp1;

    //--- extract parameters to readable variables
    if (Nparam==0) {
      d   = d_fixed;
      //--- phi is actually stored as 1-phi (confusing...)
      phi = phi_fixed;
    } else if (Nparam==1) {
      d   = param[0];
      //--- phi is actually stored as 1-phi (confusing...)
      phi = phi_fixed;
    } else {
      d   = param[0];
      //--- phi is actually stored as 1-phi (confusing...)
      phi = param[1];
    }

    //--- Allocate Memory for Hypergeometric function values
    try {
      _2F1 = new double[m];
    }
    catch (bad_alloc& ba)  {
      cerr << "bad_alloc caught: " << ba.what() << endl;
    } 

    //--- For phi=0, we have pure power-law noise
    if (fabs(phi*phi)<TINY) {
      if (d>0.75) {
        cerr << "kappa is smaller than -1.5 (" << -2.0*d << ") which is "
             << " non-stationary. Cannot use this formula." << endl;
        cerr << "1-phi: " << phi << endl;
        exit(EXIT_FAILURE);
      } 
      gamma_x[0] = tgamma(1.0-2.0*d)/pow(tgamma(1.0-d),2.0);
                                        
      for (i=0;i<m-1;i++) {
        tau = static_cast<double>(i+1);
        gamma_x[i+1] = (d+tau-1.0)*gamma_x[i]/(tau-d);
      }
 
    //--- Not pure power-law noise:
    } else { 
      //--- For d=0, _2F1 is always 1.0
      if (fabs(d)<TINY) {
        for (i=0;i<m;i++) _2F1[i] = 1.0;
      } else {
        //--- Since phi is actually stored as 1-phi, I here need to
        //    put 1- (1-phi) = phi. DONT DELETE THIS COMMENT!!!
        z       = pow(1-phi,2.0);
        k        = m-1;
        b        = d;
        a        = d   + static_cast<double>(k);
        c        = 1.0 + static_cast<double>(k);
        _2F1[m-1]= boost::math::hypergeometric_pFq( {a,b}, {c}, z);
        a       -= 1.0;
        c       -= 1.0;
        _2F1[m-2]= boost::math::hypergeometric_pFq( {a,b}, {c}, z);

        Fp1   = _2F1[m-1];
        F     = _2F1[m-2];
        for (i=m-3;i>=0;i--) {
          _2F1[i] = Fm1 = backward(a,b,c,z,F,Fp1);

          //--- prepare next round
          a  -= 1.0;
          c  -= 1.0;
          Fp1 = F;
          F   = Fm1;
        }
      }
      //--- finally, construct gamma_x
      scale = 1.0;
      for (i=0;i<m;i++) {
        gamma_x[i] = scale*_2F1[i];
        scale     *= (d+i)*(1.0-phi)/(i+1.0);
        if (std::isnan(gamma_x[i])) {
          cout << "Trouble in paradise!" << endl;
          cout << "i=" << i << ", d=" << d << ", 1-phi=" << phi << endl;
          exit(EXIT_FAILURE);
        }
      }
    } 
        
    //--- Clean up memory
    if (_2F1!=NULL) delete[] _2F1;
  }



/*! Show phi and fractional difference values
 */
//---------------------------------------------------------------------
  void GenGaussMarkov::show(double *param, double *error, double sigma)
//---------------------------------------------------------------------
  {
    double        d,T;
    Observations  &observations=Observations::getInstance();
    JSON          &json = JSON::getInstance();

    using namespace std;
    T = 1.0/(365.25*24.0*3600.0*observations.get_fs()); // T in yr

    if (std::isnan(d_fixed)) d = param[0];
    else                d = d_fixed;
    
    cout << "sigma     = " << sigma/pow(T,0.5*d)
             			<< " " << unit << "/yr^" << 0.5*d << endl;   
    
    //--- Add info to JSON file
    json.write_double("sigma",sigma/pow(T,0.5*d));

    //--- In the output write 1-phi because that's what stored 
    if (Nparam==0) {
      printf("d         = %8.4lf (fixed)\n",d_fixed);
      printf("kappa     = %8.4lf (fixed)\n",-2.0*d_fixed);
      printf("1-phi     = %le (fixed)\n",phi_fixed);
      json.write_double("d",d_fixed);
      json.write_double("kappa",-2.0*d_fixed);
      json.write_double("1-phi",phi_fixed);
    } else if (Nparam==1) {
      printf("d         = %8.4lf +/- %6.4lf\n",param[0],error[0]);
      printf("kappa     = %8.4lf +/- %6.4lf\n",-2.0*param[0],2.0*error[0]);
      printf("1-phi     = %le (fixed)\n",phi_fixed);
      json.write_double("d",param[0]);
      json.write_double("kappa",-2.0*param[0]);
      json.write_double("1-phi",phi_fixed);
    } else {
      printf("d         = %8.4lf +/- %6.4lf\n",param[0],error[0]);
      printf("kappa     = %8.4lf +/- %6.4lf\n",-2.0*param[0],2.0*error[0]);
      printf("1-phi     = %le +/- %le\n",param[1],error[1]);
      json.write_double("d",param[0]);
      json.write_double("kappa",-2.0*param[0]);
      json.write_double("1-phi",param[1]);
    }
  }



/*! Compute penalty function is parameters get out of bound
 */
//-----------------------------------------------------
  double GenGaussMarkov::compute_penalty(double *param)
//-----------------------------------------------------
  {
    double        penalty=0.0,d,phi,y;
    const double  LARGE=1.0e5,safety_factor=2.0;

    using namespace std;
    if (Nparam==0) {
      penalty = 0.0;
    } else if (Nparam==1) {
      //--- d within limits, phi is fixed
      if (param[0]>1.5) {
        penalty += (param[0]-1.5)*LARGE;
        param[0] = 1.5;
      //--- d (lower limit) If people want negative d, then I will hear it...
      } else if (param[0]<0.01) {
        penalty += (0.01-param[0])*LARGE;
        param[0] = 0.01;
      }
    } else {
      d   = param[0];
      phi = param[1];

      //cout << "IN   d=" << d << ", 1-phi=" << phi << endl;

      //--- 2F1 still explodes for large d and small 1-phi. Exclude this region
      //    using a simple line in d & log10(1-phi) space.
      //
      //    f(d) = 4*d - 11   [ f(1)=-7,  f(2)=-3
      //
      if (phi>0.0)   y = log10(param[1]);
      else           y = 9.9e99;        // will not be used

      //--- Check if log(1-phi) is below the line (2F1 is too large)
      if (y < 4.0*d - 11.0 + safety_factor) {
        penalty += ((4.0*d-11.0+safety_factor) - y)*LARGE;
        param[1] = pow(10,4.0*d-11.0+safety_factor);
      //--- param[1] is always 1-phi. The following limit is rarely encountered
      } else if (phi>0.999) {
        penalty += (phi-0.999)*LARGE;
        param[1] = 0.999;
      //--- The following limit is most critical, stay away from zero!!
      //    Another complication is that at phi=0, you have power-law and
      //    then d_max=0.5. Thus, a jump down from 2. Allowing this is 
      //    asking for trouble. I put lower limit to 1.0e-6.
      } else if (phi<1.0e-6) {
        penalty += (1.0e-6-phi)*LARGE;
        param[1] = 1.0e-6;

        //--- Now, check if 1-phi is still above the line
        y = log10(param[1]);
        if (y < 4.0*d - 11.0 + safety_factor) {
          penalty += ((4.0*d-11.0 + safety_factor) - y)*LARGE;
          param[1] = pow(10,4.0*d-11.0 + safety_factor);
        }
      }
      //--- d (upper limit)
      if (d>1.5) {
        penalty += (d-1.5)*LARGE;
        param[0] = 1.5;
      }
      //--- d (lower limit) If people want negative d, then I will hear it...
      if (d<0.01) {
        penalty += (0.01-d)*LARGE;
        param[0] = 0.01;
      }
    }
    //cout << "OUT  d=" << param[0] << ", 1-phi=" << param[1] << endl << endl;

    return penalty;
  }



/*! Set noise parameters for PSD computation.
 *
 * If phi_fixed is set, only d can vary and Nparam==1. 
 * If d_fixed is also set, nothing can vary and Nparam==0.
 * Otherwise, set both d and phi
 */
//---------------------------------------------------------------
  void GenGaussMarkov::set_noise_parameters(double *params_fixed)
//---------------------------------------------------------------
  {
    using namespace std;
    if (Nparam>=1) {
      if (std::isnan(d_fixed)==true) {
        cout << "Enter parameter value of d: ";
        cin >> d_fixed;
        params_fixed[0] = d_fixed;
      }
    }
    if (Nparam==2) {
      if (std::isnan(phi_fixed)==true) {
        cout << "Enter parameter value of 1-phi: ";
        cin >> phi_fixed;
        params_fixed[1] = phi_fixed;
      }
    }
  }



/*! Compute PSD for given frequency
 */
//-----------------------------------------------
  double GenGaussMarkov::compute_G(double lambda)
//-----------------------------------------------
  {
    return 1.0/pow(4.0*(1-phi_fixed)*pow(sin(0.5*lambda),2.0) + 
						pow(phi_fixed,2.0),d_fixed); 
  }



/*! Compute impulse response: h
 */
//---------------------------------------------------------------
  void GenGaussMarkov::compute_impulse_response(int m, double* h)
//---------------------------------------------------------------
  {
    int    i;
    double I,*params_fixed=NULL;

    //--- get values of noise parameters
    if (Nparam>0) params_fixed = new double[Nparam];
    set_noise_parameters(params_fixed);

    //--- Use Langbein (2004) formula to compute h's
    h[0] = 1.0;
    for (i=1,I=1.0;i<m;i++,I+=1.0) 
			h[i] = (d_fixed+I-1.0)/I*h[i-1]*(1.0-phi_fixed);

    if (params_fixed!=NULL) delete[] params_fixed;
  }
