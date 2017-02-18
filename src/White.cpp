/*! \file    White.cpp
 *  \author  Machiel Bos
 *
 * This extremely simple noise model is convenient for testing.
 *
 * \date 16/1/2012
 */
//==============================================================================
  #include "White.h"
  #include "Control.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <cstring>

//==============================================================================
// Subroutines
//==============================================================================


//---!!-------------
  White::White(void)
//---!!-------------
  {
    Control   &control=Control::getInstance();

    using namespace std;
    try {
      control.get_string("PhysicalUnit",unit);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    Nparam = 0;
  }



/*! The covariance matrix is a unit matrix
 */
//-----------------------------------------------------------------
  void White::get_covariance(double * /*param*/, int m, double *gamma_x)
//-----------------------------------------------------------------
  {
    int       i;

    gamma_x[0] = 1.0;
    for (i=1;i<m;i++) gamma_x[i]=0.0;
  }



/*! Nothing to show
 */
//------------------------------------------------------------
  void White::show(double * /*param*/, double * /*error*/, double sigma)
//------------------------------------------------------------
  {
    using namespace std;

    cout << "sigma     = " << sigma << " " << unit << endl;
    cout << "No noise parameters to show" << endl;
  }



/*! There is no penalty
 */
//--------------------------------------------
  double White::compute_penalty(double * /*param*/)
//--------------------------------------------
  {
    return 0.0;
  }



/*! Nothing can be set or needs to be stored for white noise
 */
//------------------------------------------------------
  void White::set_noise_parameters(double * /*params_fixed*/)
//------------------------------------------------------
  {
    // Nothing needs to be stored
  }



/*! Compute PSD for given frequency
 */
//--------------------------------------
  double White::compute_G(double /*lambda*/)
//--------------------------------------
  {
    return 1.0;
  }



/*! Compute impulse response: h
 */
//------------------------------------------------------
  void White::compute_impulse_response(int m, double* h)
//------------------------------------------------------
  {
    memset(h,0,m*sizeof(double)); //--- zero padding
    h[0] = 1.0;
  } 

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
