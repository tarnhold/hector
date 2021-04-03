/*! \file   GenGaussMarkov.cpp
 *  \author Machiel Bos
 *
 * This class implements the Generalized Gauss-Markov noise model that is
 * preferred by John Langbein Simon Williams.
 *
 * \date 5/6/2012  Santa Clara
 */
//==============================================================================
  #include "GenGaussMarkov.h"
  #include <iostream>
  #include <ostream>
  #include <gsl/gsl_sf.h>

  #define TINY 1.0e-9
//  #define DEBUG

//==============================================================================
// Subroutines
//==============================================================================


//---!!-------------------------------
  GenGaussMarkov::GenGaussMarkov(void)
//---!!-------------------------------
  {
    Control          *control=Control::getInstance();

    using namespace std;
    Nparam = 2;
    //--- Sanity check (Observations.cpp already checks if firstdifference
    //    keyword exists in the first place.
    try {
      if (control->get_bool("firstdifference")==true) {
        cerr << "Have not implemented first differenced Gauss-Markov!" << endl;
        exit(EXIT_FAILURE);
      }
    }
    catch (const char* str) {
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
    HyperGeo     Func;
    int          i,k;
    double       scale,tau,phi,a,b,c,d,z,*_2F1=NULL,Fm1,F,Fp1;

    //--- extract parameters to readable variables
    d   = param[0];
    phi = param[1];

    //--- Allocate Memory for Hypergeometric function values
    try {
      _2F1 = new double[m];
    }
    catch (bad_alloc& ba)  {
      cerr << "bad_alloc caught: " << ba.what() << endl;
    } 

    //--- For phi=0, we have pure power-law noise
    if (fabs(phi*phi)<TINY) {
      gamma_x[0] = exp(gamma(1.0-2.0*d))/pow(exp(gamma(1.0-d)),2.0);
                                        
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
        z       = pow(1-phi,2.0);
        k        = m-1;
        b        = d;
        a        = d   + static_cast<double>(k);
        c        = 1.0 + static_cast<double>(k);
        _2F1[m-1]= real(Func.compute(a,b,c,z));
        a       -= 1.0;
        c       -= 1.0;
        _2F1[m-2]= real(Func.compute(a,b,c,z));

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
        if (isnan(gamma_x[i])) {
          cout << "Trouble in paradise!" << endl;
          cout << "i=" << i << ", d=" << d << ", phi=" << phi << endl;
          exit(EXIT_FAILURE);
        }
      }

    } 
        
#ifdef DEBUG
    cout << "d=" << d << ", phi=" << phi << endl;
    for (i=0;i<m;i++) cout << "i=" << i << ", 2F1[i]=" << _2F1[i] << endl;
#endif

    //--- Clean up memory
    if (_2F1!=NULL) delete[] _2F1;
  }



/*! Show phi and fractional difference values
 */
//----------------------------------------------------
  void GenGaussMarkov::show(double *param, double *error)
//----------------------------------------------------
  {
    using namespace std;
    printf("d   = %8.4lf +/- %6.4lf\n",param[0],error[0]);
    printf("phi = %8.4lf +/- %6.4lf\n",param[1],error[1]);
  }



/*! Compute penalty function is parameters get out of bound
 */
//--------------------------------------------------
  double GenGaussMarkov::compute_penalty(double *param)
//--------------------------------------------------
  {
    double        penalty=0.0;
    const double  LARGE=1.0e3;

    using namespace std;
    //--- phi
    if (param[1]>0.9999) {
      penalty += (param[1]-0.9999)*LARGE;
      param[1] = 0.9999;
    } else if (param[1]<0.0) {
      penalty += (0.0-param[1])*LARGE;
      param[1] = 0.0;
    }
    //--- d (upper limit)
    if (fabs(param[1])>0.01) {
      if (param[0]>1.0e5) {
        penalty += (param[0]-1.0e5)*LARGE;
        param[0] = 1.0e5;
      }
    } else {
      if (param[0]>0.499) {
        penalty += (param[0]-0.499)*LARGE;
        param[0] = 0.499;
      }
    }
    //--- d (lower limit)
    if (param[0]<-0.499) {
      penalty += (-0.499-param[0])*LARGE;
      param[0] = -0.499;
    }
    //cout << "out) d=" << param[0] << " phi=" << param[1] << endl;
    return penalty;
  }



/*! Setup parameters for PSD computation
 */
//---------------------------------
  void GenGaussMarkov::setup_PSD(void)
//---------------------------------
  {
    using namespace std;
    cout << "Enter parameter value of d: ";
    cin >> d_fixed;
    cout << "Enter parameter value of phi: ";
    cin >> phi_fixed;
  }



/*! Compute PSD for given frequency
 */
//--------------------------------------------
  double GenGaussMarkov::compute_G(double lambda)
//--------------------------------------------
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
    double I;

    //--- get values of noise parameters
    setup_PSD();

    //--- Use Langbein (2004) formula to compute h's
    h[0] = 1.0;
    for (i=1,I=1.0;i<m;i++,I+=1.0) 
			h[i] = (d_fixed+I-1.0)/I*h[i-1]*(1.0-phi_fixed);
  }
