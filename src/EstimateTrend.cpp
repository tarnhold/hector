/*! \file    EstimateTrend.cpp
 *  \author  Machiel Bos
 *
 * This program estimates a trend from the observations.
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
//=============================================================================
  #include "Minimizer.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <ctime>

//  #define TIME

//=============================================================================
// Main program
//=============================================================================

  int main(int argc, char *argv[])
  {
    using namespace std;

    //--- Open correct control file
    if (argc==1) {
      Control &control = Control::getInstance("estimatetrend.ctl");
      (void) control;
    } else if (argc==2) {
      Control &control = Control::getInstance(argv[1]);
      (void) control;
    } else {
      cerr << "correct usage: estimatetrend [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Start estimatetrend
    cout << endl 
         << "************************************" << endl 
         << "    estimatetrend, version " << VERSION << "." << 
         SUBVERSION << endl << "************************************" << endl;

    Minimizer      minimizer;
#ifdef TIME
    time_t         start,end;

    time(&start);
#endif

    minimizer.solve();

#ifdef TIME
    time(&end);
    cout << "Total computing time: " << difftime(end,start) << " sec" << endl;
#endif

    //--- DesignMatrix is not a Meyer singleton so delete manually
    DesignMatrix   *designmatrix=DesignMatrix::getInstance();
    designmatrix->resetInstance();

    return EXIT_SUCCESS;
  }

/* vim:set shiftwidth=2 softtabstop=2 expandtab: */
