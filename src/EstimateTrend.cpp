/*! \file    EstimateTrend.cpp
 *  \author  Machiel Bos
 *  \version 1.0
 *
 * The main program that simply calls the minimizer class. 
 *
 * \date 16/1/2012  Santa Clara
 */
//=============================================================================
  #include "Minimizer.h"
  #include <iostream>
  #include <ostream>
  #include <cstdlib>
  #include <ctime>

//=============================================================================
// Main program
//=============================================================================

  int main(int argc, char *argv[])
  {
    Control   *control;

    using namespace std;
    //--- Open correct control file
    if (argc==1) {
      control = Control::getInstance("estimatetrend.ctl");
    } else if (argc==2) {
      control = Control::getInstance(argv[1]);
    } else {
      cerr << "correct usage: estimatetrend [controlfile.ctl]" << endl;
      exit(EXIT_FAILURE);
    }

    //--- Start estimatetrend
    cout << endl 
         << "************************************" << endl 
         << "    estimatetrend, version " << VERSION << endl    
         << "************************************" << endl;

    Minimizer      minimizer;
    time_t         start,end;

    time(&start);
    minimizer.solve();
    time(&end);
    cout << "Total computing time: " << difftime(end,start) << " sec" << endl;

    return EXIT_SUCCESS;
  }
