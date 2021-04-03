/*! \file    White.cpp
 *  \author  Machiel Bos
 *
 * This extremely simple noise model is convenient for testing.
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
  void White::get_covariance(double *param, int m, double *gamma_x)
//-----------------------------------------------------------------
  {
    int       i;

    gamma_x[0] = 1.0;
    for (i=1;i<m;i++) gamma_x[i]=0.0;
  }



/*! Nothing to show
 */
//------------------------------------------------------------
  void White::show(double *param, double *error, double sigma)
//------------------------------------------------------------
  {
    using namespace std;
    cout << "sigma     = " << sigma << " " << unit << endl;
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



/*! Nothing can be set or needs to be stored for white noise
 */
//------------------------------------------------------
  void White::set_noise_parameters(double *params_fixed)
//------------------------------------------------------
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



/*! Compute impulse response: h
 */
//------------------------------------------------------
  void White::compute_impulse_response(int m, double* h)
//------------------------------------------------------
  {
    memset(h,0,m*sizeof(double)); //--- zero padding
    h[0] = 1.0;
  } 
