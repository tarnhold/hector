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
         << "    estimatetrend, version " << VERSION << endl    
         << "************************************" << endl;

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
