/*! \file   Likelihood.cpp
 *  \author Machiel Bos
 *
 * The core class that computes the likelihood function. The exact method 
 * for computing the likelihood and Least-Squares bit is defined in the derived
 * classes AmmarGrag and FullCov classes which again make use of MLEBase.
 *
 * This class is a Meyer singleton. To combine it with the factory design
 * pattern I use the variable *method which points to the right method.
 * Nevetheless, it implies adding some wrapper subroutines. However, these are
 * easy to understand and to maintain.
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
  #include "Likelihood.h"
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


//--!!------------------------
  Likelihood::Likelihood(void)
//--!!------------------------
  {
    using namespace std;
    int           Ngaps;
    double        fraction;
    Control       &control = Control::getInstance();
    Observations  &observations = Observations::getInstance();

    Ngaps = observations.number_of_gaps();
    try {
      control.get_string("LikelihoodMethod",methodname);
    }
    //--- If keyword is not found, then use AmmarGrag is percentage of 
    //    missing data is less than 50%, otherwise use FullCov.
    catch (exception &e) {
      fraction = static_cast<double>(Ngaps)/
           static_cast<double>(observations.number_of_observations());
      if (fraction<0.5) {
        methodname = "AmmarGrag";
      } else {
        methodname = "FullCov";
      }
    }
    cout << endl << "----------------" << endl << "  " << methodname
         << endl << "----------------" << endl;
    if (methodname.compare("AmmarGrag")==0) {
      method = new AmmarGrag();
    } else if (methodname.compare("FullCov")==0) {
      method = new FullCov();
    } else {
      cerr << "Unkown Likelihood Method: " << methodname << endl;
      exit(EXIT_FAILURE);
    }
  }



/*! Free memory
 */
//---!!------------------------
  Likelihood::~Likelihood(void)
//---!!------------------------
  {
    delete method;
  }



/*! Force new initialisation of classes
 */
//-----------------------------------
  void Likelihood::reset_method(void)
//-----------------------------------
  {
    using namespace std;
    delete method;
    if (methodname.compare("AmmarGrag")==0) {
      method = new AmmarGrag();
    } else if (methodname.compare("FullCov")==0) {
      method = new FullCov();
    } else {
      cerr << "Unkown Likelihood Method: " << methodname << endl;
      exit(EXIT_FAILURE);
    }
  }
 


/*!
 */
//----------------------------------------------------
  void Likelihood::compute_LeastSquares(double *param)
//----------------------------------------------------
  {
    method->compute_LeastSquares(param);
  }



/*!
 */
//-----------------------------------------
  double Likelihood::compute(double *param)
//-----------------------------------------------------
  {
    return method->compute(param);
  }



/*!
 */
//----------------------------------------
  void Likelihood::show_leastsquares(void)
//----------------------------------------
  {
    method->show_leastsquares();
  }



/*!
 */
//-------------------------------------------------
  void Likelihood::compute_L_and_ICs(double *param)
//-------------------------------------------------
  {
    method->compute_L_and_ICs(param);
  }
