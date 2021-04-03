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

//==============================================================================
// Subroutines
//==============================================================================


//---!!-------------
  White::White(void)
//---!!-------------
  {
    Nparam = 0;
  }



/*! The covariance matrix is a unit matrix
 */
//-----------------------------------------------------------------
  void White::get_covariance(double *param, int m, double *gamma_x)
//-----------------------------------------------------------------
  {
    int       i;
    Control   *control=Control::getInstance();

    //--- Remember that covariance matrix changes when first difference is
    //    applied.
    if (control->get_bool("firstdifference")==true) {
      gamma_x[0] =  2.0;
      gamma_x[1] = -1.0;
      for (i=2;i<m;i++) gamma_x[i]=0.0;
    } else {
      gamma_x[0] = 1.0;
      for (i=1;i<m;i++) gamma_x[i]=0.0;
    }
  }



/*! Nothing to show
 */
//----------------------------------------------
  void White::show(double *param, double *error)
//----------------------------------------------
  {
    using namespace std;
    cout << "No noise parameters to show" << endl;
  }



/*! There is no penalty
 */
//--------------------------------------------
  double White::compute_penalty(double *param)
//--------------------------------------------
  {
    return 0.0;
  }



/*! Nothing needs to be stored for white noise
 */
//---------------------------
  void White::setup_PSD(void)
//---------------------------
  {
    // Nothing needs to be stored
  }



/*! Compute PSD for given frequency
 */
//--------------------------------------
  double White::compute_G(double lambda)
//--------------------------------------
  {
    return 1.0;
  }
